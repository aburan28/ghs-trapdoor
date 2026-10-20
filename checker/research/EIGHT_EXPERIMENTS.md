# Eight descent experiments

This follow-up to `DESCENT_LEADS.md` executes finite, declared experiments.
“Exhaustive” below applies to the stated finite sets, not all curves, covers,
fields, algorithms, or cryptographic families. A failed implementation or an
unresolved normalization is recorded separately from a mathematical exclusion.

All experiments use exact arithmetic. Reports include input and source hashes,
runtime versions, bounds, setup/evaluation timing, controls, and unimplemented
steps. Scalar labels may validate an evaluator but may not be its input.

| Experiment | Finite matrix and controls | Required certificate |
| --- | --- | --- |
| 1. Guided isogenies | All generic j-invariants over F343 and F2197, both quadratic twists; exclude j=0,1728; rational degree-3 first steps. Recover the existing nine-vertex component. | Actual torsion action, explicit isogenies, verified odd-prime subgroup images; count isogeny classes and Frobenius j-orbits separately. |
| 2. Trace-zero transport | Opposite-sign kernel 810; all 373 source inputs if a non-enumerative inversion is obtained; otherwise record the exact missing map and test the available trace-zero divisors. | Nonzero transfer and Frobenius anti-equivariance; abstract equal-order identification is insufficient. |
| 3. Lift signs | Matching kernel 153 and opposite kernel 810; all nonrational F49 point differences and rational point-difference controls. | Actual divisor classes, odd-order checks, plus/minus Frobenius identities; distinguish intrinsic Kummer descent from coordinates in a nonrational theta basis. |
| 4. Divisor inversion | First source/kernel 153; degree-3 divisors in a monogenic affine chart, six calibrated Weil values, fixed deterministic test inputs. | Polynomial equations, chart denominators, finite solution set, re-evaluation and sign checks; solver budgets and unsupported charts are reported. |
| 5. Obstruction filter | Both sources and every one of their twelve paired constructions. | A normalized chi18 value, including the differential factor, or an explicit unresolved result; compare with the full polynomial control without using it as the predictor. |
| 6. Dimension four | Every 2-Lagrangian for all three elliptic 2-torsion types and all eleven genus-2 branch cycle types; the existing genus-4 cover is a Jacobian control. | Kernel bases, Frobenius and pairing checks; actual curve and subgroup arithmetic for the control; no-factor screening does not certify the Jacobian locus. |
| 7. Pryms | Squarefree f(u)=u^10+a*u+b over F3, a in F3, b nonzero; C:y^2=f(x^2), B:v^2=u*f(u), D:y^2=f(u). | Free involution (x,y)->(-x,-y), genera 9/5/4, full zeta factors and odd-prime transfer controls; account for the factor two in the induced polarization. |
| 8. Wild covers | All reduced cubic classes a*x^3+b*x+c over F16 and F64, modulo h^2+h; c represents one of two constant classes. Explicit rank-two fiber-product control. | Class normalization, Frobenius spans, ramification/genus, complete quotient Frobenius polynomials, explicit cover maps and norm/conorm degree identities. |

Experiments 4 and 3 come first; experiments 1 and 5 improve screening, and
experiment 2 depends on a usable inversion and a correct twist identification.
Experiments 6–8 have separate geometry and subgroup gates.

Setup, map evaluation, search, relation generation, and linear algebra are
separate cost categories. An unimplemented category has no measured cost and
cannot be silently treated as zero. Generic square-root work is a comparison
scale, not a measured end-to-end DLP benchmark. No practical speedup is claimed.

## Results

All eight tracks have executed experiments. Several central algorithmic goals
remain unresolved; positive controls are distinguished from general algorithms.

1. **Guided search and a second reconstruction case.** The census exhausts
   5,072 generic-j models and quadratic twists: 682 over F343 and 4,390 over
   F2197. All selected negative first-step searches close, giving 558 routes
   with explicit odd-prime subgroup witnesses (78 and 480 respectively).
   These are models/routes, not independent families. See
   [census](results/eight-isogeny-search.json).

   One F2197 floor model has a prime-order 2,137 subgroup. Its degree-three
   quotient has coefficients A=1835, B=1161 in the census polynomial basis.
   All 2,137 images and generator additions are checked, and all twelve mixed
   kernels on this target are reconstructed as smooth plane quartics over F13,
   with all 28 bitangents per curve verified. Full point counts give **three
   matching and nine opposite-sign targets**, rather than the F7 example's
   eleven-to-one split. Their Jacobian orders are 2,137 and 2,259. This is a
   second characteristic and isogeny class, still at toy size; no divisor
   evaluator for the F13 targets is supplied. See
   [reconstruction](results/eight-cross-characteristic-targets.json).

2. **Trace-zero transport.** All 187 classes modulo sign of the kernel-810
   trace-zero subgroup match actual quotient-section evaluations. An adjacent
   quotient evaluation fixes the sign. **All 373 source inputs lift and all
   373 generator additions pass** in actual Jacobian arithmetic. With the
   stated generator D=[P-Frob7(P)], a consistent source-generator image is
   103D. Every image therefore satisfies Frob7(Q)=-Q. This constructs a finite
   subgroup evaluator by quotient sections plus an **O(r) inverse table**;
   it is neither scalable inversion nor a global algebraic quotient-map
   formula. See [certificate](results/eight-trace-zero-transport.json).

3. **Lift signs.** Every F49 curve point on both control quartics is tested:
   100 controls total, comprising 16 rational controls and 84 nonrational
   point differences. All 42 nonrational differences on kernel 810 have
   exact order 373. On the matching quartic the corresponding orders are
   15, 21 or 105. Frobenius fixes the rational controls and negates the point
   differences. The unordered class pair certifies intrinsic Kummer descent;
   an F7 coordinate embedding and a general normalized-map sign certificate
   are still absent. See [sign controls](results/eight-lift-signs.json).

4. **Direct inversion.** The symmetric monogenic system is implemented using
   U(x)=x^3+u2*x^2+u1*x+u0 and Z=V(x), deg(V)<=2. Common denominators, determinant
   ratios and norms supply the six equations, together with the curve
   remainder and denominator saturation. Three construction attempts exhaust
   their declared 120/300/300-second budgets; these are implementation limits,
   not nonexistence results. A second formulation introduces auxiliary values
   to avoid expanding squared determinant products. It constructs 99 equations
   in 96 variables in about one second. Known supports are used only for
   forward validation against independent Hess evaluations, never supplied to
   the solver. Every solve runs in an isolated process with a 120-second
   budget. See [sparse attempts](results/eight-divisor-inversion-circuit.json).
   No scalable inverse, global sign resolution, or larger-field scaling result
   is established.

5. **Obstruction filter.** All 24 original theta constructions are recomputed.
   The tested product features coincide: the product and its chosen square
   root have minimal field degree nine, and the product is a square there.
   This feature group contains all **22 matching and two opposite-sign**
   constructions, so these product features cannot distinguish the controls.
   The missing weight-18 differential factor can change square class under a
   twist (delta^(-27) for delta=3 over F7). This does not rule out a correctly
   normalized chi18 filter; that implementation remains open. See
   [probe](results/eight-obstruction-filter.json).

6. **Dimension four.** For each pairing model, all 200,787 RREF four-planes
   are examined and exactly 2,295 Lagrangians found. All fourteen declared
   Frobenius-module types are screened, with explicit bases for every retained
   kernel. Sage independently checks the kernels and realizes all three
   elliptic and eleven genus-2 branch types on actual curves. This exhausts
   the stated level-two module screen, not the Jacobian locus. Separately, the
   saved genus-4 cover now has an **actual nonzero order-17 divisor image**
   under D -> pi^*D+Frob7(pi^*D), with source and target orders checked. General
   Schottky filtering and quotient reconstruction remain open. See
   [kernels](results/eight-dimension-four-kernels.json) and
   [subgroup control](results/eight-dimension-four-control.json).

7. **Four-dimensional Pryms.** All six squarefree polynomials
   f(u)=u^10+a*u+b over F3, a in F3 and b nonzero, are tested. The covers
   C:y^2=f(x^2) have genus nine; B:v^2=u*f(u) has genus five and
   D:y^2=f(u) genus four. Counts through F(3^9) independently verify
   P_C=P_B*P_D in every case. The involution (x,y)->(-x,-y) is free, giving
   a four-dimensional principally polarized Prym; the restriction of the
   covering Jacobian polarization is twice principal. For f=u^10+u+2 an
   **order-41 divisor is explicitly transferred** from Jac(D) into the Prym.
   Its norm cancels as an actual divisor. Matching arbitrary restriction
   factors to Pryms and cheaper relation algorithms remain open. See
   [zeta controls](results/eight-prym.json) and
   [subgroup control](results/eight-prym-subgroup.json).

8. **Artin--Schreier classes and the classical baseline.** The F16/F64 sweep
   exhausts **17,408 reduced presentations** across the two declared families,
   with both constant trace classes. It records constant-field obstructions,
   Frobenius spans and conductor genera. These are classes with a fixed
   x-coordinate; x-scaling identifies further ordinary presentations.
   Normalizing the two-pole family gives 438 comparisons over all proper
   subfields: **382 regular cases have exactly the classical checker's genus**;
   56 irregular cases remain uncertified as descents over the smaller field.
   No lower genus than this regular baseline is found. Sage independently
   verifies 1,014 genus-one counts and 4,352 field multiplication pairs.

   The rank-two cubic control has the F2 model z^4+z=x^3. Complete genus-three
   counts verify its three elliptic zeta factors over F16. A degree-two
   quotient to y^2+y=a*x^3 over F4 gives an actual nonzero order-three pullback
   and F2-rational trace divisor. This control is supersingular. General
   ordinary norm/conorm descent maps remain open. See
   [class sweep](results/eight-wild-classes.json),
   [classical comparison](results/eight-wild-classical-comparison.json), and
   [subgroup control](results/eight-wild-subgroup.json).

## Reproduction

From `checker/`, install the pinned reconstruction backend. Native commands:

The complete pipeline is also available as
`python -m research.run_eight_experiments`; `--list` prints its exact command
arrays and `--sage-command` selects a Sage command prefix. Each stage archives
its command, log, return code and wall time. `--manifest-only` hashes an already
completed result set.

```bash
python -m research.wild_artin_schreier --output research/results/eight-wild-classes.json
python -m research.wild_classical_comparison --output research/results/eight-wild-classical-comparison.json
python -m research.dimension_four_kernels --output research/results/eight-dimension-four-kernels.json
python -m research.obstruction_filter_probe --output research/results/eight-obstruction-filter.json
python -m research.quotient_reconstruction --export-divisor-inputs --kernel-index 810 --output research/results/eight-trace-zero-inputs.json
python -m research.eight_cost_baseline --output research/results/eight-cost-baseline.json
```

The following require Sage 10.6. Run the isogeny census before the native
cross-characteristic reconstruction, and run the independent checks last:

```bash
sage -python research/eight_geometry_sage.py --experiment 1 --output research/results/eight-isogeny-search.json
python -m research.cross_characteristic_targets --output research/results/eight-cross-characteristic-targets.json
sage -python research/eight_geometry_sage.py --experiment 3 --output research/results/eight-lift-signs.json
sage -python research/eight_geometry_sage.py --experiment 6 --output research/results/eight-dimension-four-control.json
sage -python research/eight_geometry_sage.py --experiment 7 --output research/results/eight-prym.json
sage -python research/eight_geometry_sage.py --experiment 8 --output research/results/eight-wild-cover.json
sage -python research/cover_subgroups_sage.py --kind prym --output research/results/eight-prym-subgroup.json
sage -python research/cover_subgroups_sage.py --kind wild --output research/results/eight-wild-subgroup.json
sage -python research/trace_zero_transport_sage.py --output research/results/eight-trace-zero-transport.json
sage -python research/divisor_inversion_system_sage.py --budget 300 --output research/results/eight-divisor-inversion-det3.json
sage -python research/divisor_inversion_circuit_sage.py --budget 120 --output research/results/eight-divisor-inversion-circuit.json
sage -python research/eight_independent_checks_sage.py --output research/results/eight-independent-checks.json
python -m unittest discover -v
```

Sage is separate from Python CI. The official container used here is
`sagemath/sagemath:10.6`, digest
`sha256:19995db6194f4a4bab18ce9a88556fd15b9ed5e916b4504fefe618a7796ddbdb`.
Earlier inversion-construction source snapshots are archived with their
reports. Timeouts are wall-clock engineering bounds, not formal complexity
lower bounds. Research commands reject optimized Python that removes asserts.

The final local suite has **71 passing Python tests**. Changed Python files
pass Ruff and the diff passes whitespace checks. A repository-wide Ruff run
also reports 28 pre-existing findings, all verified to lie in files unchanged
from the base commit; those unrelated files were not rewritten. Sage checks
are recorded separately, including all 36 extension-field counts for the twelve
F13 targets, all 2,137 isogeny additions, the fourteen dimension-four module
realizations/checks, and the Artin–Schreier field/count oracles. Remote CI has
not run for this local branch.

## Cost and evidence limits

The [toy BSGS baseline](results/eight-cost-baseline.json) solves the specified
order-373 elliptic DLP using 20 baby entries and 26 counted group additions.
Field-table setup is timed separately. Different arithmetic backends and
representations prevent treating these toy timings as a fair cryptographic
benchmark. The trace-zero inverse uses 187 actual Jacobian classes modulo
sign and remains O(r); the sparse systems have not replaced it.

Reports distinguish setup, evaluation and total elapsed time where measured.
Unmeasured phases and unimplemented relation generation/linear algebra are
not zero-cost phases. No complete attack-cost comparison, practical speedup,
new vulnerable cryptographic family, or changed NIST result is established.
