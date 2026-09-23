#!/usr/bin/env python
"""Run the packaged parallel Tailer module from a source checkout."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from Tailer.tailer_parallel import main


if __name__ == "__main__":
    main()
