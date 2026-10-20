# GHS research tracking

Last updated: 2026-10-20.

| Work | Pull request | Evidence and current status |
| --- | --- | --- |
| Classical binary GHS checker and bounded explicit degree-3 scans | [#2](https://github.com/aburan28/ghs-trapdoor/pull/2) | Draft; 16 local tests pass. Exact classical parameters, published presets and explicit isogeny transport; no complete descent/DLP implementation. |
| Published genus-2 Richelot audit | [#3](https://github.com/aburan28/ghs-trapdoor/pull/3) | Draft, stacked on #2. Three BK first-step obstructions and 45 generic/Kummer quotient constructions excluded for the selected weak-cover branch condition. |
| Explicit genus-2 pair and genus-4 cover | [#3](https://github.com/aburan28/ghs-trapdoor/pull/3) | Toy construction; polynomial map identity and full Frobenius polynomials verified. Combined suite: 25 local tests pass. No DLP or practical speedup claimed. |
| Longer genus-2 isogeny paths | Follow-up to #3 | Untested. Requires Frobenius-stable kernels even when individual points/factors are not rational. |
| Odd-degree genus-2 isogenies | Follow-up to #3 | Untested. Absence of rational Richelot kernels does not exclude these. |
| Other target spaces and field connections | [#4](https://github.com/aburan28/ghs-trapdoor/pull/4) | Draft; primary-source catalogue of established methods and explicit untested directions. |
| Rational quotients of a cubic elliptic Weil restriction | [#4](https://github.com/aburan28/ghs-trapdoor/pull/4) | Kernel screening implemented at ell=2/3. Published example: zero rational (2,2,2) kernels; 16 rational 3-torsion Lagrangians, 12 passing the no-coordinate-factor test. Combined suite: 34 local tests pass. |
| Explicit quotient reconstruction | [#6](https://github.com/aburan28/ghs-trapdoor/pull/6), stacked on #5 | All 24 candidate kernels give smooth plane quartics over F7. 22 match the required Frobenius polynomial; two opposite-sign targets have order 315, obstructing an F7-rational order-373 map. |
| Explicit quotient-to-Kummer maps | [#6](https://github.com/aburan28/ghs-trapdoor/pull/6) | All 373 source points evaluated per curve; exactly 187 images, with only sign pairs colliding. All 81 tested kernel translations preserve the image for each source. |
| Jacobian divisor lift | [#6](https://github.com/aburan28/ghs-trapdoor/pull/6) | Completed Sage experiment: all 373 points, 373 generator additions and the scalar-37 relation checked per source model. The second uses verified Galois transport into the first target. Divisor inversion remains O(r). |
| Algebraic Kummer-to-divisor inversion component | [#8](https://github.com/aburan28/ghs-trapdoor/pull/8), stacked on #7 | Eight-function diagonal system: rank 112, kernel 8. Three degree-12 fibers solved, each giving two effective divisors; exact residual identities and independent FLINT checks. No Jacobian-group enumeration. Calibration from elliptic quotient coordinates and exceptional charts remain pending. |
| Conditional Frobenius lift-sign certificates | [#8](https://github.com/aburan28/ghs-trapdoor/pull/8) | Twelve actual base-field Jacobian divisor tests distinguish the two known twist possibilities: eleven positive signs, one negative at kernel 810. No extension-field curve counts used to form the certificates; not a general equal-order isogeny test or χ18 implementation. |
| CM-guided candidates in additional isogeny classes | [#8](https://github.com/aburan28/ghs-trapdoor/pull/8) | Six Frobenius j-orbits over F13^3 and F19^3, in three elliptic isogeny classes. All independently recounted, prime-order generators checked, twelve mixed 3-kernel candidates each. Four first-class isogeny neighbors lose the condition. New quotient targets remain unconstructed. |
| Wild Artin–Schreier extension classes | [#8](https://github.com/aburan28/ghs-trapdoor/pull/8) | Sixteen equivalent degree-12 presentations reduced to a genus-2 model over F8; 8192 affine point transports checked. Base Jacobian order 65 cannot preserve order-37/109 extension subgroups. Supersingular toy control; multiple-pole covers and actual norm-map tests remain pending. |
| Isogeny routes and endomorphism orders | [#7](https://github.com/aburan28/ghs-trapdoor/pull/7), stacked on #6 | Exact nine-vertex F343 component: three vertices have 12 mixed 3-kernel candidates, six have zero and reach the earlier example via one 3-isogeny plus conjugation. All 6,714 edge point images and generator additions checked; Sage independently verifies all edges and the two CM class polynomials. |
| Opposite-sign target over a quadratic extension | [#7](https://github.com/aburan28/ghs-trapdoor/pull/7); [new notes](research/DESCENT_LEADS.md) | Actual order-373 divisor in ker(Frob7+1) on the kernel-810 quartic Jacobian over F49. No explicit elliptic-source map into that target or relation-generation algorithm is constructed. Its class modulo sign is F7-rational, explaining why Kummer rationality alone is insufficient. |
| Broader experiments and lower-cost filters | [#7](https://github.com/aburan28/ghs-trapdoor/pull/7) | Six prioritized directions include CM-guided parameter generation, Frobenius lift signs, nonenumerative inversion, χ18, dimension-four Jacobian/Prym tests and Artin–Schreier extension classes. Untested steps are labeled explicitly. |
| Actual cubic-field curve census | [#5](https://github.com/aburan28/ghs-trapdoor/pull/5), stacked on #4 | 16,140 nonsingular equations checked. Two sampled prime-order F343 inputs, with j outside F7, have zero rational 2-kernels and 12 mixed 3-kernel candidates each. Quotient targets are not reconstructed. |
| Explicit Prym factors | [#5](https://github.com/aburan28/ghs-trapdoor/pull/5) | 32 genus-3 covers with elliptic/genus-2 quotients; full Frobenius polynomial factorization verified in every case. 23 cases have an odd prime in the genus-2 factor order absent from the elliptic factor order. |
| Torus DLP transfer controls | [#5](https://github.com/aburan28/ghs-trapdoor/pull/5) | Four nodal curves: all 440 group-law pairs and 40 scalar logs checked. The same candidate map fails on all four smooth elliptic perturbations. |

Related existing implementation: [#1 — Sage trap-door demo and Lean formalization](https://github.com/aburan28/ghs-trapdoor/pull/1). Its reported results were inspected for context; Sage/Lean validation was not rerun in this session.

Review order: #2 targets `main`; #3 targets `codex/classical-ghs-checker`. After #2 is merged, retarget #3 to `main` and check its diff. #4 targets `codex/genus2-weil-descent-audit`; after #3 is merged, retarget #4 to `main` and check its diff. #5 targets `codex/alternative-descent-targets`; after #4 is merged, retarget #5 to `main` and check its diff. #6 targets `codex/descent-experiment-batch`; after #5 is merged, retarget #6 to `main` and check its diff. #7 targets `codex/explicit-quotient-reconstruction`; after #6 is merged, retarget #7 to `main` and check its diff. #8 targets `codex/isogeny-descent-leads`; after #7 is merged, retarget #8 to `main` and check its diff. PRs #2–#8 are drafts; no merge has been performed.

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

The [new batch notes](research/EXPERIMENT_BATCH.md) explain the search bounds, actual prime-order curve inputs, full Prym example, and positive/negative torus controls. At the #5 stage, the local suite had 45 tests, including independent finite-field comparisons and group-law checks of both Frobenius eigenspaces via quadratic twists. The batch also fixes derivative trimming in the squarefreeness helper when the characteristic divides the polynomial degree.

No previously unreported cryptographic vulnerability or NIST binary-curve break has been established.

The [quotient reconstruction notes](research/QUOTIENT_RECONSTRUCTION.md) give equations, exact twist obstructions and the map methodology. The optional reconstruction backend is pinned in `research/requirements-reconstruction.txt`; the #6-stage suite had 55 passing Python tests when it was installed.

A source-equivalence audit finds one underlying Weil-restriction case: the two inputs are related by coefficient conjugation and an explicit scaling. All 373 point transports and the twelve-kernel correspondence are verified. Counts above refer to constructions, not independent curve families.

Completed divisor-map evidence: [Sage report](research/results/quotient-divisor-verification.json), with input SHA256, explicit divisor witnesses, shared-target kernel correspondence and all subgroup checks.

The #7 [descent leads](research/DESCENT_LEADS.md) brought the suite to **59 passing Python tests**. Four new tests cover the graph, independent kernel-point-sum oracles in characteristics 5,7,11, invalid kernels and trace-zero evidence boundaries. The separate Sage 10.6 report records 6,714 independently evaluated edge images, class-polynomial checks, and a trace-zero divisor witness; no new practical speedup is established.

The [next follow-up](research/DESCENT_FOLLOWUP.md) brings the full suite to **67 passing Python tests**. Three new Sage 10.6 experiments completed: CM selection, conditional sign certificates and symmetric Kummer inversion. Reports and tests distinguish the new inversion component from a complete nonenumerative source-to-Jacobian transfer, which is not yet implemented. No change to the NIST conclusions follows.
