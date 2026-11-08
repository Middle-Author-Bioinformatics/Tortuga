#!/usr/bin/env python3
"""Run ITSME, MitoLooper, or both using one pipeline selector."""
import sys
from mito_pipeline import main

if __name__=="__main__":
    try:
        sys.exit(main(require_pipeline=True))
    except (ValueError,RuntimeError,OSError,KeyError) as error:
        print(f"ERROR: {error}",file=sys.stderr)
        sys.exit(1)
