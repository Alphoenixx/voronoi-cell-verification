# Machine verification for *Asymptotically Optimal Online Selection of a Large Voronoi Cell*

`verify.py` is the script accompanying the manuscript

> A. Sadhu, *Asymptotically Optimal Online Selection of a Large Voronoi Cell*,

which extends S. Har-Peled, *The Prophet and the Voronoi Diagram*, ESA 2026,
LIPIcs vol. 388, Article 41, [`doi:10.4230/LIPIcs.ESA.2026.41`](https://doi.org/10.4230/LIPIcs.ESA.2026.41).

It re-derives every computational claim in the paper from scratch, and is written
so that a referee can see at a glance which claims it **proves** and which it
merely **tests**. That distinction is the point of the script, and it is enforced
in the output labels, not just in this file.

## Running it

```bash
pip install sympy mpmath numpy scipy
python -u verify.py
```

Roughly 3 minutes. The script uses half the logical cores; override with
`VERIFY_PROCS=n`. Output is flushed line by line, so it can be watched as it
runs. Exit status is non-zero if any check fails.

Current run: **64 checks, 0 failures**, about 2½ minutes on 8 processes.

## What is proved, and what is only tested

**Proof-grade.** Exact rational or algebraic arithmetic, or verified `mpmath`
interval arithmetic. No floating-point quantity is load-bearing anywhere in this
group.

| Claim | How it is settled |
|---|---|
| Lemma 4.4 | the closed form for `f`, cross-checked against the raw definition of `D_p` by ray shooting, using no formula from the paper |
| Thms 4.5 / 4.6 / 4.7 | all eleven inequalities, decided exactly in `sympy` |
| Cor 4.3 | `f(1/2,1/2) = 4√2/3 − 5/3`, exactly |
| Prop 4.9 | both `f(q)` and `∇f(q)` enclosed in verified interval arithmetic at all 64 probes, so the supporting halfplane is valid for the *true* gradient; outward rounding to rationals; polygon intersection, the cyclic **order** of the vertices (sign tests on `Fraction`s, not `atan2`, then certified convex-CCW), and the shoelace area, all in exact rationals |
| Remark 4.8 | the fundamental-domain identity, in exact rationals |
| Lemma 5.4 | the one-hit and two-hit areas, in closed form |
| Thm 5.6 | the two-cone expansion `1 − 2(1−a)^s + (1−2a)^s = s(s−1)a² + O(s³a³)` behind the universal reliability barrier — including that the linear term vanishes, which is why a single hit is not enough |
| Section 6 | the area-transfer counting step of Lemma 6.7 and the arithmetic `r log r ≤ (1+o(1)) log m` behind the architecture bound of Prop 6.8; `∫₀^{2π} g = 1` and `1/8 ≤ g ≤ 1/4` for the polar area density `g(θ) = 1/(8max(cos²θ, sin²θ))` of the tile, in `sympy`; the exact `8/r` and `16/r` neighbour separations; the satellite containment in `(1+2/r)Q` for `r ≥ 4`; and the polynomial identity `r(r−8)² − (r−20)(r+2)² = 140r + 80`, which gives the gap `r(1−8/r)²/(1+2/r)² ≥ r − 20` for every `r > 0` |

**Consistency tests.** Monte Carlo or simulation. These are evidence, not proof,
and are labelled as such in the output.

- Lemma 4.1 — concavity of `f^(1/2)` on random segments
- Lemmas 2.1 / 2.2 — the flower inequality, and that a Voronoi flower is empty
- Remark 4.10 — the two numerical constants
- Lemma 5.1 — the scaling of `E[Y]`
- Thm 6.1 — the equal-area hidden-tile construction is built independently at `r = 21, 41, 81` and its `r − O(1)` area gap measured; the measured deficits are `4.79, 5.24, 5.47`, bounded and well inside the proved constant `20`
- Prop 2.5 / Thm 3.2 — simulated torus point sets

On that last row: `2^d n A / log n → 1` carries a correction of order
`log log n / log n`, so at any simulable `n` the measured value is near 2, not
near 1. The script prints the ratio
`(2^d n A/log n − 1) / (log log n/log n)`, which is ≈ 4 and near-constant across
`n = 8000, 20000, 40000` — so the observed value is the predicted one. It is
checked for magnitude and trend, never for equality.

**Not machine-checked.** The probabilistic arguments themselves: the conditioning
in the lower bound of Theorem 5.5, and the posterior on the unresolved gadgets in
Theorem 6.1. Those are hand proofs. `ALL CHECKS PASSED` does not cover them, and
the script says so.

## License

MIT. See [LICENSE](LICENSE).
