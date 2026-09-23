#!/usr/bin/env python
"""Parallel Tailer: 3' end tail analysis with a multiprocessing pool.

Reimplements the *global* (GTF annotation + SAM/BAM) mode of
TimNicholsonShaw/tailer using ``multiprocessing`` so a single input file can use
many cores. Output CSV columns and per-row values are identical to upstream
Tailer, so downstream parsing does not need to change.

Why a rewrite instead of threads: upstream Tailer is pure-Python and GIL-bound;
the hot loop is per-read gffutils/SQLite gene-overlap queries. Threads would not
help. This version keeps every alignment of a read together and collapses reads
by sequence *before* distributing work, so the multimapper "best alignment"
logic is preserved exactly while the gene queries run across processes.

Usage (drop-in for global mode):
    python tailer_parallel.py -a ref/gencode.v49.annotation.gtf -p 16 \
        rawfiles/WT_4SU_uniq.bam

    # multiple files (each parallelized internally, processed one at a time)
    python tailer_parallel.py -a ref/hg38-tRNAs.withnames.gtf -p 16 \
        rawfiles/*_uniq.bam

Writes <input_basename>_tail.csv next to each input, matching upstream Tailer.
"""

import argparse
import csv
import os
import sqlite3
import multiprocessing as mp
import time

import pysam
import gffutils
from Bio.Seq import Seq

try:
    from tqdm import tqdm
except ImportError:  # tqdm is an upstream dep but keep this runnable without it
    def tqdm(x, **_):
        return x


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# --- Picklable alignment record ------------------------------------------------
# pysam AlignedSegment objects cannot cross process boundaries cheaply, so we
# extract only the fields Tailer needs into a plain tuple.
class AlnRec:
    __slots__ = (
        "query_sequence",
        "query_length",
        "is_reverse",
        "reference_name",
        "pos",
        "reference_length",
        "cigartuples",
    )

    def __init__(self, aln):
        self.query_sequence = aln.query_sequence
        self.query_length = aln.query_length
        self.is_reverse = aln.is_reverse
        self.reference_name = aln.reference_name
        self.pos = aln.pos
        self.reference_length = aln.reference_length
        self.cigartuples = aln.cigartuples

    # Enable pickling of __slots__ objects for multiprocessing.
    def __getstate__(self):
        return {s: getattr(self, s) for s in self.__slots__}

    def __setstate__(self, state):
        for k, v in state.items():
            setattr(self, k, v)


def getOrMakeGTFdb(GTForDB):
    """Build or load a gffutils SQLite db of gene features (one time).

    Mirrors upstream Tailer: reduces the GTF to gene rows, builds <pre>.db, and
    reuses it on subsequent runs. Returns the db path so each worker can open
    its own connection.
    """
    pre, ext = os.path.splitext(GTForDB)
    if ext == ".db":
        return GTForDB
    if os.path.exists(pre + ".db"):
        return pre + ".db"

    with open(GTForDB, "r") as gtffile, open(pre + "_temp.gtf", "w") as tempfile:
        for line in gtffile:
            if line.startswith("#"):
                continue
            if line.split("\t")[2] == "gene":
                tempfile.write(line)

    print("Creating GTF database...")
    gffutils.create_db(
        pre + "_temp.gtf",
        dbfn=pre + ".db",
        force=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )
    os.remove(pre + "_temp.gtf")
    return pre + ".db"


def getHandleOnBam(samOrBamFile):
    pre, ext = os.path.splitext(samOrBamFile)
    if ext.lower() in (".fasta", ".fastq"):
        raise Exception("FASTA/FASTQ can only be used in local mode.")
    if ext.lower() == ".sam":
        return pysam.AlignmentFile(samOrBamFile, "r")
    elif ext.lower() == ".bam":
        return pysam.AlignmentFile(samOrBamFile, "rb")
    raise Exception(ext + " is an unsupported file type")


def collapseBySequence(handledBAM, read=2):
    """Group alignments by read name, then collapse reads by sequence.

    Returns a dict: sequence -> [count, [AlnRec, ...]]. All alignments of a read
    stay together (needed for correct best-alignment selection); reads that
    share an identical sequence are pooled and their count incremented, exactly
    as upstream Tailer does.
    """
    name_to_alns = {}
    for aln in handledBAM:
        if read == 2 and aln.is_read1:
            continue
        if read == 1 and aln.is_read2:
            continue
        if aln.is_unmapped:
            continue
        name_to_alns.setdefault(aln.query_name, []).append(AlnRec(aln))

    seq_dict = {}
    for _, alns in name_to_alns.items():
        seq = alns[0].query_sequence
        entry = seq_dict.get(seq)
        if entry is None:
            seq_dict[seq] = [1, alns]
        else:
            entry[0] += 1  # duplicate sequence: just bump the count
    name_to_alns.clear()  # drop per-read query_name strings; seq_dict shares the AlnRecs
    return seq_dict


# --- Worker globals (set once per process via Pool initializer) ----------------
_DB = None
_REVCOMP = False
_THRESH = 100
_SEQOUT = False
_MATURE = 0
_SE_READ1 = False  # single-end: interpret reads in read1 (True) or read2 (False) orientation


def _open_immutable_db(db_path):
    """Open the gffutils DB read-only with SQLite immutable mode.

    immutable=1 disables locking and journaling, so many worker processes can
    read the same prebuilt .db concurrently without 'database is locked' errors
    on shared/NFS filesystems (e.g. Biowulf /vf/users).
    """
    uri = "file:" + os.path.abspath(db_path) + "?immutable=1"
    conn = sqlite3.connect(uri, uri=True, check_same_thread=False)
    try:
        return gffutils.FeatureDB(conn)
    except Exception:
        # Older gffutils that only accept a path: build, then swap in the
        # immutable connection for the actual (heavy, concurrent) queries.
        db = gffutils.FeatureDB(db_path)
        db.conn = conn
        return db


def _init_worker(db_path, rev_comp, thresh, seq_out, mature, se_read1):
    global _DB, _REVCOMP, _THRESH, _SEQOUT, _MATURE, _SE_READ1
    _DB = _open_immutable_db(db_path)
    _REVCOMP = rev_comp
    _THRESH = thresh
    _SEQOUT = seq_out
    _MATURE = mature
    _SE_READ1 = se_read1


def _eff_reverse(is_reverse):
    """Effective read orientation after the flip toggles.

    The gene-region strand lookup and the 3'-end branch must agree; a flip from
    --rev_comp or single-end read1 has to move BOTH, or the lookup queries the
    opposite strand and finds no overlapping gene (dropping the read).
    """
    return is_reverse != (_REVCOMP != _SE_READ1)


def _make_tail(a, gene):
    """Replicates upstream Tailer.Tail: (geneID, gene_name, threeEnd, tailLen, tailSeq)."""
    try:
        gene_name = gene["gene_name"][0]
    except KeyError:
        gene_name = None
    try:
        geneID = gene["gene_id"][0]
    except KeyError:
        geneID = None

    # Effective orientation carries the --rev_comp / single-end read1 flip.
    positive = _eff_reverse(a.is_reverse)

    if positive:
        read_end = a.pos + a.reference_length
        threeEnd = read_end - gene.stop  # pysam 0-indexed, GTF 1-indexed
        if a.cigartuples[-1][0] == 4:  # 3' softclip
            tailLen = a.cigartuples[-1][1]
            tailSeq = a.query_sequence[a.query_length - tailLen:]
        else:
            tailLen = 0
            tailSeq = None
    else:
        threeEnd = gene.start - a.pos - 1
        if a.cigartuples[0][0] == 4:  # 5' softclip, needs rev-comp
            tailLen = a.cigartuples[0][1]
            tailSeq = reverse_complement(a.query_sequence[:tailLen])
        else:
            tailLen = 0
            tailSeq = None

    return (geneID, gene_name, threeEnd, tailLen, tailSeq, _eff_reverse(a.is_reverse))


def _row_for_seq(seq, count, alns):
    tails = []
    for a in alns:
        # Query the strand matching the read's effective orientation (post-flip).
        strand = "+" if _eff_reverse(a.is_reverse) else "-"
        for gene in _DB.region(
            seqid=a.reference_name,
            start=a.pos,
            end=a.pos + a.reference_length,
            strand=strand,
        ):
            tails.append(_make_tail(a, gene))

    if not tails:
        return None

    best = min(abs(t[2]) for t in tails)
    bestTails = [t for t in tails if abs(t[2]) == best]

    ensIDs = "|".join(sorted(str(t[0]) for t in bestTails))
    geneNames = "|".join(sorted(str(t[1]) for t in bestTails))

    b0 = bestTails[0]
    threeEnd, tailLen, tailSeq, is_reverse = b0[2], b0[3], b0[4], b0[5]

    if threeEnd < -_THRESH or threeEnd > _THRESH:
        return None

    end_pos = threeEnd + tailLen - _MATURE

    if not _SEQOUT:
        return [count, ensIDs, geneNames, end_pos, tailLen, tailSeq]

    out_seq = seq if is_reverse else reverse_complement(seq)
    return [out_seq, count, ensIDs, geneNames, end_pos, tailLen, tailSeq]


def _process_chunk(chunk):
    rows = []
    for seq, count, alns in chunk:
        row = _row_for_seq(seq, count, alns)
        if row is not None:
            rows.append(row)
    return rows


def _iter_chunks(seq_dict, size):
    """Yield chunks of (seq, count, alns), popping from seq_dict as we go.

    Draining the source dict lets the parent release each read's AlnRecs once its
    chunk has been dispatched, instead of holding the whole file plus separate
    work/chunks copies in memory at once.
    """
    chunk = []
    while seq_dict:
        seq, (cnt, alns) = seq_dict.popitem()
        chunk.append((seq, cnt, alns))
        if len(chunk) >= size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _process_chunks_with_recovery(chunk_iter, nchunks, processes, db_path,
                                  rev_comp, thresh, seq_out, mature, se_read1,
                                  chunk_timeout):
    """Process chunks in a bounded pool and recover stalled in-flight chunks."""
    all_rows = []
    pending = {}
    indexed_chunks = enumerate(chunk_iter)
    exhausted = False
    max_inflight = max(1, processes * 2)

    pool = mp.Pool(
        processes=processes,
        initializer=_init_worker,
        initargs=(db_path, rev_comp, thresh, seq_out, mature, se_read1),
        maxtasksperchild=4,
    )

    def submit_until_full():
        nonlocal exhausted
        while not exhausted and len(pending) < max_inflight:
            try:
                idx, chunk = next(indexed_chunks)
            except StopIteration:
                exhausted = True
                break
            pending[idx] = (chunk, pool.apply_async(_process_chunk, (chunk,)))

    submit_until_full()
    completed = 0
    last_progress = time.monotonic()

    try:
        with tqdm(total=nchunks) as progress:
            while pending:
                made_progress = False
                for idx, (chunk, result) in list(pending.items()):
                    if not result.ready():
                        continue
                    try:
                        rows = result.get()
                    except Exception as exc:
                        raise RuntimeError(
                            "Worker failed while processing chunk {}".format(idx)
                        ) from exc
                    all_rows.extend(rows)
                    del pending[idx]
                    completed += 1
                    progress.update()
                    made_progress = True
                    last_progress = time.monotonic()
                    submit_until_full()

                if made_progress:
                    continue

                if time.monotonic() - last_progress > chunk_timeout:
                    stalled = sorted(pending)
                    print(
                        "No worker result returned for {} seconds after {}/{} "
                        "chunks. Recovering {} in-flight chunks in the parent: {}"
                        .format(chunk_timeout, completed, nchunks, len(stalled), stalled),
                        flush=True,
                    )
                    pool.terminate()
                    pool.join()
                    _init_worker(db_path, rev_comp, thresh, seq_out, mature, se_read1)
                    for idx in stalled:
                        chunk, _ = pending[idx]
                        all_rows.extend(_process_chunk(chunk))
                        progress.update()
                    for _, chunk in indexed_chunks:
                        all_rows.extend(_process_chunk(chunk))
                        progress.update()
                    return all_rows

                time.sleep(1)
    except Exception:
        pool.terminate()
        pool.join()
        raise

    print("Collected all worker results; stopping worker pool.", flush=True)
    pool.terminate()
    pool.join()
    return all_rows


def runFile(bam_file, db_path, processes, rev_comp, thresh, seq_out, mature, se_read1,
            chunk_timeout):
    pre, _ = os.path.splitext(bam_file)
    out_loc = pre + "_tail.csv"

    print("Tailing file: " + bam_file)
    with getHandleOnBam(bam_file) as handle:
        seq_dict = collapseBySequence(handle)

    n = len(seq_dict)
    # Small chunks reduce per-worker memory and make tail-end stragglers less costly.
    size = max(1, n // (max(1, processes) * 16))
    nchunks = (n + size - 1) // size if n else 0
    print("Calculating tails for {} unique sequences on {} processes...".format(
        n, processes))

    all_rows = []
    if processes <= 1:
        _init_worker(db_path, rev_comp, thresh, seq_out, mature, se_read1)
        for chunk in tqdm(_iter_chunks(seq_dict, size), total=nchunks):
            all_rows.extend(_process_chunk(chunk))
    else:
        all_rows = _process_chunks_with_recovery(
            _iter_chunks(seq_dict, size),
            nchunks,
            processes,
            db_path,
            rev_comp,
            thresh,
            seq_out,
            mature,
            se_read1,
            chunk_timeout,
        )
    del seq_dict

    key_col = 1 if seq_out else 0  # sort by Count
    print("Sorting {} output rows...".format(len(all_rows)), flush=True)
    all_rows.sort(key=lambda x: x[key_col], reverse=True)

    if seq_out:
        header = ["Sequence", "Count", "EnsID", "Gene_Name",
                  "End_Position", "Tail_Length", "Tail_Sequence"]
    else:
        header = ["Count", "EnsID", "Gene_Name",
                  "End_Position", "Tail_Length", "Tail_Sequence"]

    print("Writing {} output rows to {}...".format(len(all_rows), out_loc), flush=True)
    with open(out_loc, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for row in all_rows:
            writer.writerow(row)

    print("Wrote " + out_loc + " to disk.")


def main():
    parser = argparse.ArgumentParser(
        description="Parallel Tailer: SAM/BAM to tail file using a GTF (global mode)")
    parser.add_argument("-a", "--annotation", required=True, metavar="",
                        help="GTF annotation (or prebuilt gffutils .db)")
    parser.add_argument("files", nargs="+",
                        help="SAM or BAM formatted files")
    parser.add_argument("-p", "--processes", type=int, default=os.cpu_count(),
                        metavar="", help="Number of worker processes (default: all cores)")
    parser.add_argument("-t", "--threshold", type=int, default=100, metavar="",
                        help="Max distance from mature end to include (default=100)")
    parser.add_argument("-r", "--rev_comp", action="store_true",
                        help="Reverse complement reads")
    parser.add_argument("-s", "--sequence", action="store_true",
                        help="Output nucleotide sequences to file")
    parser.add_argument("-read", "--read", type=int, default=2, metavar="",
                        help="Paired-end only: which read holds the 3' end (default=2)")
    parser.add_argument("-se", "--se_read", type=int, default=2, choices=(1, 2), metavar="",
                        help="Single-end only: interpret reads in read1 or read2 orientation "
                             "(default=2, matching paired-end -read default; use 1 for "
                             "PEAR-merged reads that carry the R1 orientation)")
    parser.add_argument("-m", "--mature", type=int, default=0, metavar="",
                        help="Mature end adjustment")
    parser.add_argument("--chunk-timeout", type=int, default=900, metavar="",
                        help="Seconds to wait for the next worker chunk before failing "
                             "instead of hanging forever (default=900)")
    args = parser.parse_args()

    db_path = getOrMakeGTFdb(args.annotation)

    for bam_file in args.files:
        runFile(
            bam_file,
            db_path,
            processes=max(1, args.processes),
            rev_comp=args.rev_comp,
            thresh=args.threshold,
            seq_out=args.sequence,
            mature=args.mature,
            se_read1=(args.se_read == 1),
            chunk_timeout=args.chunk_timeout,
        )


if __name__ == "__main__":
    main()
