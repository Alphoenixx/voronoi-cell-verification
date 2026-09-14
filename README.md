# Machine verification for *Asymptotically Optimal Online Selection of a Large Voronoi Cell*

This repository accompanies

> A. Sadhu, *Asymptotically Optimal Online Selection of a Large Voronoi Cell*,

which extends S. Har-Peled, *The Prophet and the Voronoi Diagram*, ESA 2026,
LIPIcs vol. 388, Article 41, [`doi:10.4230/LIPIcs.ESA.2026.41`](https://doi.org/10.4230/LIPIcs.ESA.2026.41).

`verify.py` is the computational core. `verify_submission.py` is the recommended
submission-facing runner: it executes the same core without changing any
mathematical computation and translates legacy draft labels in the printed
output to the numbering of the final manuscript. `MANUSCRIPT_NUMBERING.md`
records the complete label map.

## Running it

```bash
pip install sympy mpmath numpy scipy
python -u verify_submission.py
```

Running `python -u verify.py --no-pause` performs the identical checks but may
print result numbers from the earlier draft.

The suite takes roughly 3 minutes on a typical multicore machine. It uses half
the logical cores by default; override this with `VERIFY_PROCS=n`. Output is
flushed line by line, and the exit status is non-zero if any check fails.

Current reference run: **65 checks, 0 failures**, about 2½ minutes on 8
processes.

## What is proved, and what is only tested

**Proof-grade.** Exact rational or algebraic arithmetic, or verified `mpmath`
interval arithmetic. No floating-point approximation is used in a proof-critical
comparison.

| Claim | How it is settled |
|---|---|
| Lemma 5.4 | the closed form for `f`, cross-checked against the raw definition of `D_p` by ray shooting, using no formula from the paper |
| Thms 5.5 / 5.6 / 5.7 | all eleven inequalities, decided exactly in `sympy` |
| Cor. 5.3 | `f(1/2,1/2) = 4√2/3 − 5/3`, exactly |
| Prop. 5.9 | both `f(q)` and `∇f(q)` enclosed in verified interval arithmetic at all 64 probes, so the supporting halfplane is valid for the *true* gradient; outward rounding to rationals; polygon intersection, the cyclic **order** of the vertices (sign tests on `Fraction`s, not `atan2`, then certified convex-CCW), and the shoelace area, all in exact rationals |
| Remark 5.8 | the fundamental-domain identity, in exact rationals |
| Lemma 6.4 | the extremal case: a bisector through the centre of `B` leaves exactly half, in `sympy`. The lemma itself is a symmetry argument in the paper and needs no segment area; the script's segment and slab fractions are illustration, listed below |
| Thm. 6.6 | the two-cone expansion `1 − 2(1−a)^s + (1−2a)^s = s(s−1)a² + O(s³a³)` behind the universal reliability barrier, including the vanishing linear term |
| Section 7 | the area-transfer counting step of Lemma 7.10 and the arithmetic `r log r ≤ (1+o(1)) log m` behind the architecture bound of Prop. 7.11; `∫₀^{2π} g = 1` and `1/8 ≤ g ≤ 1/4` for the square's polar area density; the exact neighbour-separation calculations; the satellite containment estimate; and the polynomial identity `r(r−8)² − (r−20)(r+2)² = 140r + 80` |

**Consistency tests.** Monte Carlo or simulation. These are evidence, not proof,
and are labelled as such in the output.

- Lemma 5.1 — concavity of `f^(1/2)` on random segments
- Lemmas 3.1 / 3.2 — the flower inequality and the empty-flower property
- Remark 5.10 — the two numerical constants
- Lemma 6.1 — the scaling of `E[Y]`
- Thm. 7.4 — the equal-area hidden-tile construction is built independently at `r = 21, 41, 81` and its `r − O(1)` area gap measured; the measured deficits are `4.79, 5.24, 5.47`, bounded and well inside the proved constant `20`
- Prop. 3.5 / Thm. 4.2 — simulated torus point sets

On that last row, `2^d n A / log n → 1` carries a correction of order
`log log n / log n`, so at simulable values of `n` the measured value is near 2,
not near 1. The script prints
`(2^d n A/log n − 1) / (log log n/log n)`, which is approximately 4 and
near-constant across `n = 8000, 20000, 40000`. This is checked for magnitude and
trend, never for equality.

**Not machine-checked.** The probabilistic arguments themselves: the conditioning
in the lower bound of Theorem 6.5, the generalized isolation tradeoff of
Theorem 4.2, and the geometric/posterior proof of the hidden-tile Theorem 7.4.
Those are hand proofs. `ALL CHECKS PASSED` does not cover them, and the script
says so.

## License

MIT. See [LICENSE](LICENSE).
