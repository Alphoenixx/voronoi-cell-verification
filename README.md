# Machine verification for *Asymptotically Optimal Online Selection of a Large Voronoi Cell*

This repository accompanies

> A. Sadhu, *Asymptotically Optimal Online Selection of a Large Voronoi Cell*,

which extends S. Har-Peled, *The Prophet and the Voronoi Diagram*, ESA 2026,
LIPIcs vol. 388, Article 41, DOI 10.4230/LIPIcs.ESA.2026.41.

The repository is synchronized directly with the numbering of the final
submission manuscript. There is no legacy-number translation layer.

`verify.py` contains the computational checks. `verify_submission.py` is the
submission-facing launcher named in the manuscript and simply executes
`verify.py`; it performs no relabelling.

## Running it

```bash
pip install sympy mpmath numpy scipy
python -u verify_submission.py --no-pause
```

Equivalently:

```bash
python -u verify.py --no-pause
```

The suite uses half the logical cores by default; override this with
`VERIFY_PROCS=n`. Output is flushed line by line and the exit status is non-zero
if any check fails.

Current reference run: **65 checks, 0 failures**, about 2.5 minutes on 8
processes.

## What is proved, and what is only tested

### Proof-grade checks

These use exact rational/algebraic arithmetic or verified `mpmath` interval
arithmetic. No floating-point approximation is used in a proof-critical
comparison.

| Manuscript result | What the script certifies |
|---|---|
| Theorems 5.5, 5.6, 5.7 | The eleven square-cell inequalities, exactly in `sympy` |
| Corollary 5.3 | `f(1/2,1/2) = 4*sqrt(2)/3 - 5/3`, exactly |
| Proposition 5.9 | Verified enclosures of `f(q)` and `grad f(q)` at all 64 probes; outward rational rounding; exact polygon intersection, cyclic order, and shoelace area |
| Remark 5.8 | Fundamental-domain area identity, exactly |
| Lemma 6.4 | The extremal one-hit bisector-through-the-centre case leaves exactly one half |
| Theorem 6.6 | The two-cone expansion in (6.5), including the vanishing linear term |
| Lemma 7.10 | The area-transfer counting step in (7.13) |
| Proposition 7.11 | The arithmetic behind (7.14) and the random-arrival averaging bound |
| Section 7 | Polar-density normalization, neighbour-separation arithmetic, satellite-containment estimate, and the exact gap polynomial used for (7.5)--(7.8) |

### Consistency tests

These are numerical cross-checks or simulations and are not used as proofs.

- Lemma 5.4: the closed form for `f` is cross-checked against the raw definition
  of `D_p` by numerical ray shooting.
- Lemma 5.1: concavity of `f^(1/2)` on random segments.
- Lemmas 3.1 and 3.2: flower inequality and empty-flower property on sampled
  instances.
- Remark 5.10: the two numerical constants.
- Lemma 6.1: scaling of `E[Y]`.
- Lemma 6.4: illustrative segment/slab fractions beyond the exact extremal check.
- Theorem 7.4: the equal-area construction is built independently for
  `r = 21, 41, 81` and its `r - O(1)` area gap is measured.
- Proposition 3.5 / Theorem 4.2: simulated torus point sets.

For the last row, `2^d n A / log n -> 1` has a correction of order
`log log n / log n`, so at simulable values of `n` the measured value can be near
2 rather than near 1. The script therefore checks the magnitude and trend of the
second-order normalized quantity, not equality.

### Not machine-checked

The probabilistic arguments themselves are hand proofs: in particular the
conditioning in the lower bound of Theorem 6.5, the generalized isolation
tradeoff of Theorem 4.2, and the geometric/posterior proof of Theorem 7.4.
`ALL CHECKS PASSED` does not claim to certify those hand arguments.

## Repository contents

- `verify.py` -- verification suite with final-manuscript numbering.
- `verify_submission.py` -- no-translation submission launcher.
- `LICENSE` -- MIT license.

## License

MIT. See `LICENSE`.
