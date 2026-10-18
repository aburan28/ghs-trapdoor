# GHS research tracking

Last updated: 2026-10-18.

| Work | Pull request | Evidence and current status |
| --- | --- | --- |
| Classical binary GHS checker and bounded explicit degree-3 scans | [#2](https://github.com/aburan28/ghs-trapdoor/pull/2) | Draft; 16 local tests pass. Exact classical parameters, published presets and explicit isogeny transport; no complete descent/DLP implementation. |
| Published genus-2 Richelot audit | [#3](https://github.com/aburan28/ghs-trapdoor/pull/3) | Draft, stacked on #2. Three BK first-step obstructions and 45 generic/Kummer quotient constructions excluded for the selected weak-cover branch condition. |
| Explicit genus-2 pair and genus-4 cover | [#3](https://github.com/aburan28/ghs-trapdoor/pull/3) | Toy construction; polynomial map identity and full Frobenius polynomials verified. Combined suite: 25 local tests pass. No DLP or practical speedup claimed. |
| Longer genus-2 isogeny paths | Follow-up to #3 | Untested. Requires Frobenius-stable kernels even when individual points/factors are not rational. |
| Odd-degree genus-2 isogenies | Follow-up to #3 | Untested. Absence of rational Richelot kernels does not exclude these. |
| Other target spaces and field connections | [#4](https://github.com/aburan28/ghs-trapdoor/pull/4) | Draft; primary-source catalogue of established methods and explicit untested directions. |
| Rational quotients of a cubic elliptic Weil restriction | [#4](https://github.com/aburan28/ghs-trapdoor/pull/4) | Kernel screening implemented at ell=2/3. Published example: zero rational (2,2,2) kernels; 16 rational 3-torsion Lagrangians, 12 passing the no-coordinate-factor test. Combined suite: 34 local tests pass. |
| Additional target reconstruction and point transfer | Follow-up to #4 | Untested. Candidate kernels still require a suitable Jacobian over the base field, explicit maps and an attack-cost comparison. |
| Actual cubic-field curve census | [#5](https://github.com/aburan28/ghs-trapdoor/pull/5), stacked on #4 | 16,140 nonsingular equations checked. Two sampled prime-order F343 inputs, with j outside F7, have zero rational 2-kernels and 12 mixed 3-kernel candidates each. Quotient targets are not reconstructed. |
| Explicit Prym factors | [#5](https://github.com/aburan28/ghs-trapdoor/pull/5) | 32 genus-3 covers with elliptic/genus-2 quotients; full Frobenius polynomial factorization verified in every case. 23 cases have an odd prime in the genus-2 factor order absent from the elliptic factor order. |
| Torus DLP transfer controls | [#5](https://github.com/aburan28/ghs-trapdoor/pull/5) | Four nodal curves: all 440 group-law pairs and 40 scalar logs checked. The same candidate map fails on all four smooth elliptic perturbations. |

Related existing implementation: [#1 — Sage trap-door demo and Lean formalization](https://github.com/aburan28/ghs-trapdoor/pull/1). Its reported results were inspected for context; Sage/Lean validation was not rerun in this session.

Review order: #2 targets `main`; #3 targets `codex/classical-ghs-checker`. After #2 is merged, retarget #3 to `main` and check its diff. #4 targets `codex/genus2-weil-descent-audit`; after #3 is merged, retarget #4 to `main` and check its diff. #5 targets `codex/alternative-descent-targets`; after #4 is merged, retarget #5 to `main` and check its diff. PRs #2–#5 are drafts; no merge has been performed.

Reproduce from the repository root:

```bash
cd checker
python -m unittest discover -v
python research/genus2_audit.py --output /tmp/genus2-audit.json
python -m research.genus2_cover_example --output /tmp/genus2-cover.json
python -m research.weil_quotient_kernels --output /tmp/weil-quotient-kernels.json
python -m research.cubic_curve_census --output /tmp/cubic-curve-census.json
python -m research.prym_cover_experiments --output /tmp/prym-cover-experiments.json
python -m research.torus_transfer_experiments --output /tmp/torus-transfer-experiments.json
```

Results and equations: [research notes](research/README.md), [published-parameter audit](research/results/genus2-audit.json), [toy cover](research/results/toy-cover.json). Alternative targets: [research map](research/ALTERNATIVE_TARGETS.md) and [kernel results](research/results/weil-quotient-kernels.json). The Python 3.10/3.13 workflow reports remote CI status in each PR; local test results above do not imply that remote CI has completed.

Further experiments should be tracked in their own PR with exact inputs, search bounds, reproducible results, subgroup/map checks, comparison with known attacks and primary-source references. Separate demonstrated results from conjectures and untested leads.

The [new batch notes](research/EXPERIMENT_BATCH.md) explain the search bounds, actual prime-order curve inputs, full Prym example, and positive/negative torus controls. The complete local suite now has 45 tests, including independent finite-field comparisons and group-law checks of both Frobenius eigenspaces via quadratic twists. The batch also fixes derivative trimming in the squarefreeness helper when the characteristic divides the polynomial degree.

No previously unreported cryptographic vulnerability or NIST binary-curve break has been established.
