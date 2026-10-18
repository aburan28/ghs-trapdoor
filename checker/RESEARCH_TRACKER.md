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

Related existing implementation: [#1 — Sage trap-door demo and Lean formalization](https://github.com/aburan28/ghs-trapdoor/pull/1). Its reported results were inspected for context; Sage/Lean validation was not rerun in this session.

Review order: #2 targets `main`; #3 targets `codex/classical-ghs-checker`. After #2 is merged, retarget #3 to `main` and check its diff. #4 targets `codex/genus2-weil-descent-audit`; after #3 is merged, retarget #4 to `main` and check its diff. PRs #2–#4 are drafts; no merge has been performed.

Reproduce from the repository root:

```bash
cd checker
python -m unittest discover -v
python research/genus2_audit.py --output /tmp/genus2-audit.json
python -m research.genus2_cover_example --output /tmp/genus2-cover.json
python -m research.weil_quotient_kernels --output /tmp/weil-quotient-kernels.json
```

Results and equations: [research notes](research/README.md), [published-parameter audit](research/results/genus2-audit.json), [toy cover](research/results/toy-cover.json). Alternative targets: [research map](research/ALTERNATIVE_TARGETS.md) and [kernel results](research/results/weil-quotient-kernels.json). The Python 3.10/3.13 workflow reports remote CI status in each PR; local test results above do not imply that remote CI has completed.

Further experiments should be tracked in their own PR with exact inputs, search bounds, reproducible results, subgroup/map checks, comparison with known attacks and primary-source references. Separate demonstrated results from conjectures and untested leads.

No previously unreported cryptographic vulnerability or NIST binary-curve break has been established.
