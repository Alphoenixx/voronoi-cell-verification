#!/usr/bin/env python3
"""Submission-facing launcher for the final-manuscript verification suite.

All theorem, lemma, section, and equation labels live directly in verify.py.
This file performs no translation or relabelling.
"""
from pathlib import Path
import runpy

runpy.run_path(str(Path(__file__).with_name("verify.py")), run_name="__main__")
