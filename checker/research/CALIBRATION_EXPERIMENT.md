# Experiment: connect the elliptic quotient to the algebraic inverse

Status: **designed; the source-side preflight has run. Target calibration, blind inversion and seed-free initialization have not run.** This protocol tests the missing link identified in [the follow-up](DESCENT_FOLLOWUP.md). The methodology starts from [Tian, §3.2](https://arxiv.org/html/2012.07173v1#S3.SS2); success would implement that link, rather than establish a new attack.

## Question and hypotheses

Can the existing 64 quotient sections be projected to eight independent coordinates, calibrated to the verified symmetric Kummer basis, and inverted to the correct Jacobian image of a supplied elliptic point without a point-label lookup?

The main hypothesis is that, after matching the curve charts, divisor origin, characteristic labels and scalar normalizations, there is a single nonsingular 8×8 matrix M such that

    [ζ(D(P))] = [M·q(P)],

where q is the eight-coordinate quotient vector and [D(P)−B] is the actual image φ(P). The equality is projective: the scalar multiplying a vector may vary with P.

A second hypothesis is that the inverse over the required fields still has the expected generic degree-12 fiber: two effective degree-3 divisors, each with six orderings. Their classes relative to B are negatives of each other. A third is that an evaluation at P+G resolves this sign consistently using the known image of G.

Separate two outcomes:

| Profile | Initialization | What a pass establishes |
| --- | --- | --- |
| A: calibrated pilot | One certified source/image anchor inherited from the enumerative toy experiment | Correct inversion per challenge after bounded calibration; initialization still depends on the previous experiment. |
| B: seed-free initialization | Normalized target Weil functions and identity values; no inherited source/image anchor | A complete nonenumerative initialization and point-transfer path for this case, subject to the stated identity and coverage checks. |

Passing A does not imply B. Neither establishes a practical attack-speed improvement.

## Fixed inputs and current preflight

Use E: y²=x³+248x+42 over F343=F7[z]/(z³+z²+1), with encoding a0+7a1+49a2, G=(0,176), and prime order r=373. Use kernel 153, with its actual basis from the [preflight manifest](results/calibration-preflight.json). Use the corresponding first quartic from the reconstruction report.

The exact input hashes, field modulus, section vectors and challenge coordinates are frozen in that manifest. Its SHA256 is

    cfca312b80fe7a35da6cd8d31ee44b195885aad5fbd6b1de07b1fc9dfe33dc65

The [reproducible preflight](calibration_preflight.py) evaluates the 64 source sections at nG for n=1,…,32. It checks rank eight and selects section indices **0,1,2,3,4,5,6,7**. The frozen projective frame is

    G, 2G, 3G, 4G, 5G, 6G, 7G, 8G, 11G.

Every eight of these nine vectors have rank eight. The target frame has not yet been evaluated. Rank on a finite pilot is a necessary check, not a proof of a global section identity.

The manifest also freezes **32 challenge points**, sampled as coordinate pairs from the tag `ghs-calibration-holdout-v1`. Their discrete-log labels and target images were not computed. Both signs of every pilot point are excluded. Each challenge P and its auxiliary P+G also occupy distinct sign classes from every other challenge or auxiliary input. In this odd-order Weierstrass group, excluding reused x-coordinates excludes reused ± classes. The stream used 99 counters within a fixed bound of 10,000.

Do not select replacement challenge points after seeing target results. A changed protocol gets a new version and a new report.

## Gate 1: align the charts, origins and fields

There are three curve presentations. Write F for the original quartic and use column coordinates.

    Ccert: F(Mcert·v)=0,  Mcert=[[1,0,0],[0,1,1],[0,0,1]],
    Csym:  F(Msym ·v)=0,  Msym =[[5,1,1],[0,1,0],[1,0,0]].

The certificate-to-symmetric coordinate map is v_sym=Msym^(-1)·Mcert·v_cert. Verify this polynomial identity and the three anchor-point memberships before evaluating divisors. Equality of curve invariants is insufficient for this step.

In Ccert, the certified anchor is

    D*=[Q+2P0−3O],
    Q=(1:2:4), P0=(1:1:0), O=(1:4:1).

In Csym it becomes

    D*=[(4:6:3)+2(0:1:0)−3(1:5:5)].

Use this divisor itself as φ(G) for profile A. The previous scalar label `110` is a validation result and must not be used to label challenge points. An arbitrary generator of an equal-order Jacobian is not an acceptable replacement for D*.

Let H be the canonical hyperplane divisor, and D9,D19,D26 the contact divisors with the existing characteristic labels. Set

    δ=D9+D19+D26−H,  B=δ+P∞,  P∞=(1:0:0) on Csym.

Verify 2δ~K_C, hence 2B~K_C+2P∞. For a target class X, construct an effective degree-3 divisor D_X~X+B by Riemann–Roch. Evaluate ζ on D_X and interpret every recovered divisor as [D_X−B]. Do not substitute the degree-zero anchor's support directly into the three-point evaluator: its displayed representative has a negative part and repeated/infinite support.

The field handling is a real dependency. The source sections use **L=F7^18**, whereas the tested inverse uses F49 and splits its fibers in F7^12. Since 18 does not divide 12, those embeddings cannot be reused. Parameterize the inverse over L and, initially, use **Ω=F7^108**, which contains every degree-1,2,3 residue field over L. Store explicit irreducible moduli and the embedding L→Ω; verify the defining relation and Frobenius compatibility exactly. This extension contains all roots only after the degree-at-most-three divisor decomposition has been established; it is not a universal splitting-field claim for arbitrary failed systems.

An optional reduction to smaller fields must first prove descent of B, the relevant projective vectors and the equations. Frobenius-fixed coordinates may then be converted through a checked embedding. Do not replace coordinates by traces, norms or truncated coefficient lists.

## Gate 2: fit the nine-point calibration and freeze it

For the nine fixed training labels n, form X_n=nD* by Jacobian addition and choose D_n~X_n+B. Compute z_n=ζ(D_n), and use the corresponding frozen q_n. This requires nine target evaluations and does not enumerate the target group.

Check that both collections are projective frames. Solve

    (M q_n)_j (z_n)_k − (M q_n)_k (z_n)_j = 0

using a nonzero pivot k of each z_n and the other seven coordinates j. The resulting 63×64 system must have rank 63. Its one-dimensional nullspace must give det(M)≠0. Normalize M by its first nonzero entry, and store its exact coefficients and SHA256 before opening target validation data.

**Nine-point fitting is not a validation result.** Any two projective frames in P7 admit a unique projective transformation. Even a wrong origin or unrelated map can fit those nine pairs. All consequential evidence comes from unused points, divisor identities and, ultimately, the section-identity check.

If a frozen training point is exceptional, implement its declared chart treatment or record a failed feasibility gate. Do not quietly replace it with a point chosen using validation outcomes.

## Gate 3: test the fixed map on 32 unused points

Use separate interfaces:

- The reference verifier may use the old toy result to compute the expected image. Its enumeration, scalar labels and expected divisors stay outside the evaluator's input.
- The evaluator receives only public arithmetic data, M, B, the certified anchor permitted by profile A, and coordinate inputs P. Its parser rejects validation tables, discrete-log labels and expected-output fields.

First compare [M q(P)] with [ζ(D_ref(P))] for all 32 challenges by exact cross-products. Require **32/32 agreement**. This gate isolates calibration errors from errors in polynomial inversion. Do not refit M after a mismatch.

Next, on the first eight predeclared challenges, run the inverse over the checked fields. Continue to all 32 only if those eight pass exactly. The remaining 24 were already frozen and cannot be replaced. For a regular fiber require:

1. A zero-dimensional ideal of degree 12 and twelve distinct ordered solutions over the verified splitting field.
2. Exactly two unordered degree-3 divisors, each accounting for six permutations.
3. Exact curve membership and all eight projective-coordinate identities, checked independently with FLINT.
4. The residual relation D+Dbar~K_C+2P∞, certified by the conic/polynomial identity or an equally exact divisor calculation.
5. Classes [D−B] and [Dbar−B] are negatives; the recovered class matches the reference after sign resolution.

Represent closed points as divisors over their field of definition. Do not require all three support points to lie in L, or confuse a nonrational support point with a nonrational divisor class.

## Gate 4: resolve signs and test the group map

Invert the mapped vectors for P and P+G. They supply the candidate sets {X,−X} and {Y,−Y}. Select the unique X whose X+D* lies in the second set. Check the addition on actual Jacobian classes. Apart from identity degeneracies, both candidates cannot satisfy this relation in the odd-order subgroup: that would force 2X=0 or 2D*=0.

Keep P+G as a coordinate input obtained by elliptic addition. The evaluator must not use a source scalar label to choose the sign. Handle O, G and −G explicitly and record those cases separately from the blind score.

For exceptional fibers, allow only the predeclared translations P+aG for a=0,1,2. Resolve the shifted image using its successor, then subtract aD*. This can avoid a bad effective-divisor chart without changing M. Log every failed attempt. Permit at most three attempts, hence at most six inverse systems per point; caching a shared successor is allowed and must be logged. Exhausted attempts are failures of the implemented coverage, not permission to drop a challenge.

After all 32 primary points pass, run eight scalar checks φ(37P_i)=37φ(P_i), i=0,…,7, and eight addition checks

    φ(P_i+P_(i+16))=φ(P_i)+φ(P_(i+16)),  i=0,…,7.

Finally verify [373]X=0 and Frob7(X)=X for every accepted class on the matching target. Frobenius is tested on the class, not each support coordinate.

## Negative and boundary controls

| Control | Required observation |
| --- | --- |
| Swap two input coordinates after M is frozen | The altered pipeline must not pass the entire exact holdout suite. |
| Use the old origin in place of B | Exact divisor/reference comparisons must expose the untranslated constant shift. |
| Change the kernel or input-report hash | Reject inconsistent calibration metadata before evaluation. |
| Kernel-810 opposite-sign quartic | Retain the existing rationality obstruction: it cannot receive an order-373 subgroup over F7. Do not reuse M or D* on this different target. |
| Arbitrary order-373 trace-zero generator over F49 | Do not accept it as the image of G without a certified quotient correspondence. Equal order alone supplies no geometric calibration. |
| Repeated support, a point at infinity, or an all-zero raw ζ evaluation | Report a chart boundary and test the declared translations. Never classify it as a nonexistent Jacobian image. |

The opposite-sign target's actual quotient map remains a separate experiment. The known trace-zero witness does not supply its missing source/image anchor.

## Profile B: remove dependence on the old enumerative map

Only start B after A identifies and verifies the coordinate conventions. Build normalized Jacobian 2-Weil functions from the labeled bitangents and δ. Match their cocycle and labels to the quotient functions, and normalize at the Jacobian identity. Compute the ξ-to-ζ change of basis on generic target triples, where both sides can be evaluated directly. This avoids using source/image pairs or knowing φ(G).

Validate that change of basis on unused target triples. To claim a global calibration, prove the relevant section identities after clearing the common poles in the coordinate ring; agreement on the finite order-373 group alone is insufficient. This is the direct normalized-function route underlying the published method.

Then invert q(G) algebraically to obtain {D_G,−D_G}. Choose one orientation once, and resolve all other signs using that same chosen D_G. An unavailable generic anchor may be replaced only by cG for c∈{1,2,3}, with the known inverse of c modulo 373 used to recover an oriented generator image. No reference-table search is allowed. The verifier may compare with D* up to one fixed global sign; it cannot choose a different sign for each point.

Run B with a public-input manifest that excludes the old divisor-verification report, validation-point table, scalar `110`, and any stored source/image lookup. Audit file reads and record the absence of enumeration and discrete-log calls. A pass requires the same frozen 32-point and group-law checks, plus the section-identity certificate. This is the criterion for removing the enumerative initialization claim.

## Bounds, measurements and interpretation

Use the pinned FLINT backend and Sage 10.6. Freeze limits before running: 2 GiB memory, 60 CPU seconds per inverse system, and 3600 wall seconds for profile A. A timeout is an inconclusive result under that engineering budget. Report exactly which gates ran; do not turn a timeout into a mathematical impossibility claim.

Record setup, source evaluation, matrix multiplication, Gröbner basis, root solving, divisor construction, Frobenius checks and sign resolution separately. Record field degrees, equation/fiber degrees, systems attempted, cache use, exceptional cases and peak memory. Profile A must disclose the inherited anchor's enumerative origin. Validation-only oracle costs are recorded separately and cannot be hidden inside the evaluator.

Compare with the old method including its O(r) precomputation and storage. Its cached per-point lookup can be fast; comparing only online lookup time or only one side's setup cost is misleading. Report amortized totals for one point, 32 points and the full tiny subgroup, without interpreting two toy fields as an asymptotic speedup.

After B passes, reconstruct a target for the independent F13³ candidate A=563,B=436, r=2281 and repeat with a fresh calibration. That target has not yet been constructed. Its success would provide an implementation generality check; a practical DLP claim would additionally require relation generation, linear algebra and a complete comparison with existing attacks.

## Outputs and stop rules

The [machine-readable protocol](examples/calibration-experiment.json) fixes gates and budgets. A run should produce a calibration certificate, evaluator transcript, isolated validation report, controls report and cost profile, with hashes connecting them.

Primary success for A is 32/32 exact images within the fixed translation budget, all sixteen group-law checks, all class/Frobenius checks, and passing controls. B additionally needs seed-free initialization and the global section certificate. Any failed gate stops the dependent phase. Classify the failure as normalization/origin, field embedding, calibration mismatch, unsupported chart, algebraic solver, sign selection, or resource limit. Preserve the input; never weaken the criterion to equality of orders or Kummer values alone.

To reproduce only the completed preflight:

```bash
python -m research.calibration_preflight --output /tmp/calibration-preflight.json
```

That command does not run the target experiment. The previous 67-test suite covers the existing components; no new calibration success is claimed by this design.
