#!/usr/bin/env python3
"""Run the verification suite with labels synchronized to the final manuscript.

The numerical and symbolic computations live in verify.py and are not modified
here.  This runner streams its output while replacing labels from the earlier
draft by the corresponding labels in the submission version.
"""

from __future__ import annotations

import os
import subprocess
import sys


# Use placeholders so that, for example, 4.1 -> 5.1 is not subsequently turned
# into 6.1 by another replacement.
REPLACEMENTS = [
    ("Lemmas 2.1/2.2", "Lemmas 3.1/3.2"),
    ("Lemmas 2.1 and 2.2", "Lemmas 3.1 and 3.2"),
    ("Lemma 2.1", "Lemma 3.1"),
    ("Lemma 2.2", "Lemma 3.2"),
    ("Prop 2.5", "Prop 3.5"),
    ("Proposition 2.5", "Proposition 3.5"),
    ("Thm 3.2", "Thm 4.2"),
    ("Theorem 3.2", "Theorem 4.2"),
    ("Lemma 4.1", "Lemma 5.1"),
    ("Cor 4.3", "Cor 5.3"),
    ("Corollary 4.3", "Corollary 5.3"),
    ("Lemma 4.4", "Lemma 5.4"),
    ("Thm 4.5/4.6/4.7", "Thm 5.5/5.6/5.7"),
    ("Theorem 4.5", "Theorem 5.5"),
    ("Theorem 4.6", "Theorem 5.6"),
    ("Theorem 4.7", "Theorem 5.7"),
    ("Thms 4.5 / 4.6 / 4.7", "Thms 5.5 / 5.6 / 5.7"),
    ("Remark 4.8", "Remark 5.8"),
    ("Prop 4.9", "Prop 5.9"),
    ("Proposition 4.9", "Proposition 5.9"),
    ("Remark 4.10", "Remark 5.10"),
    ("Lemma 5.1", "Lemma 6.1"),
    ("Lemma 5.4", "Lemma 6.4"),
    ("Theorem 5.5", "Theorem 6.5"),
    ("Thm 5.5", "Thm 6.5"),
    ("Theorem 5.6", "Theorem 6.6"),
    ("Thm 5.6", "Thm 6.6"),
    ("Thm 6.4", "Thm 7.4"),
    ("Theorem 6.4", "Theorem 7.4"),
    ("Lemma 6.10", "Lemma 7.10"),
    ("Prop 6.11", "Prop 7.11"),
    ("Proposition 6.11", "Proposition 7.11"),
    ("Section 6", "Section 7"),
    ("formula (4.1)", "formula (5.1)"),
    ("(6.13)", "(7.13)"),
    ("(6.14)", "(7.14)"),
    ("No floating-point value is load-bearing.",
     "No floating-point approximation is used in a proof-critical comparison."),
]


def translate(line: str) -> str:
    placeholders: list[tuple[str, str]] = []
    out = line
    for i, (old, new) in enumerate(REPLACEMENTS):
        token = f"@@LABEL_{i}@@"
        if old in out:
            out = out.replace(old, token)
            placeholders.append((token, new))
    for token, new in placeholders:
        out = out.replace(token, new)
    return out


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    core = os.path.join(here, "verify.py")
    cmd = [sys.executable, "-u", core, "--no-pause", *sys.argv[1:]]
    proc = subprocess.Popen(
        cmd,
        cwd=here,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(translate(line))
        sys.stdout.flush()
    return proc.wait()


if __name__ == "__main__":
    raise SystemExit(main())
