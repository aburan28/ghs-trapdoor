# Pursuing the descent leads

This batch completes four bounded experiments: an algebraic Kummer-to-divisor inversion component, conditional Frobenius-sign certificates, CM-guided candidates in three additional elliptic isogeny classes, and Artin–Schreier presentation controls. They advance the implementation and search pipeline. They do not establish a new vulnerable cryptographic family or an attack-cost improvement.

## Algebraic inversion without enumerating the Jacobian

The earlier point transfer inverted quotient coordinates by enumerating all 187 Jacobian classes modulo sign. The new [Sage implementation](symmetric_kummer_inversion_sage.py) implements a different part of the published approach: inversion of the eight-dimensional symmetric Kummer model in [Tian, §3.2](https://arxiv.org/html/2012.07173v1#S3.SS2).

For the first reconstructed quartic, the rational bitangent is X+6Y+2Z=0. A coordinate change puts its two contacts at P∞=(1:0:0), Q∞=(0:1:0), with the bitangent Z=0. The transformed equation is

    2X²Y²−X³Z−X²YZ−XY²Z+Y³Z−2X²Z²−2Y²Z²−2Z⁴=0.

The eight functions

    1, x, y, x², xy, y², x³, x²y

form the required Riemann–Roch basis for L(2P∞+2K_C). Their pole bounds follow from the two contacts of Z=0. They are independent, as also witnessed by evaluation, and Riemann–Roch gives dimension eight for this degree-10 divisor on a genus-3 curve.

There are 120 symmetric triple products. Evaluating on pairs (P,P,Q) gives a matrix with **rank 112 and kernel dimension eight**. The implementation then proves all eight resulting functions vanish on the diagonal by exact polynomial reduction modulo the two curve equations. This supplies symbolic identities beyond the finite evaluation sample. All coefficients descend to F7.

For each of three deterministic test inputs, the eight projective values yield a zero-dimensional polynomial system of **degree 12**. A Gröbner basis followed by linear elimination and univariate roots of degree at most three recovers all twelve ordered solutions. They form exactly **two unordered effective degree-3 divisors**, with six permutations each. No scalar labels or Jacobian-group enumeration are used by this inversion.

An independent python-flint check evaluates the eight functions on every recovered divisor and checks the curve equations. It also verifies the residual relation

    D+Dbar ~ K_C+2P∞.

Specifically, a conic in the span {1,x,y,x²} cuts D+Dbar. Substituting its equation for y into the quartic gives exactly the product of the six x-root factors, including multiplicities. This is an explicit divisor-equivalence certificate. A repeated-point diagonal is a negative control: all eight values vanish there, so it is outside the supported projective chart.

**Remaining boundary:** the new solver takes coordinates in this symmetric Kummer basis. The elliptic quotient currently exports seven selected theta-section coordinates, in a different presentation. A complete transfer still needs eight independent quotient sections, origin/theta-characteristic alignment, and the linear calibration to this basis. Repeated supports and points at infinity also need additional charts. The earlier complete EC-to-Jacobian transfer remains enumerative; this batch does not silently replace its cost claim.

The [saved inversion report](results/symmetric-kummer-inversion-sage.json) contains the coordinate change, basis matrix, test inputs, twelve-solution fiber counts, recovered supports and residual conics. The small precomputation enumerates points of C(F49); it never enumerates the Jacobian. The fixed splitting field F7^12 contains the degree-at-most-three residue fields over F49 needed by these tests. Results are limited to these three fibers, with a bounded triangular solver rather than a general variety solver.

## Conditional sign certificates before extension-field point counts

For the twelve kernels on the first source model, the geometric quotient construction leaves two possible base-field Jacobian classes, with polynomials T⁶±29T³+343 and group orders 373 or 315. The two models from the original census are Galois-related, so the second is not counted again.

The new [Sage sign test](quotient_sign_certificates_sage.py) chooses distinct F7-points P,O and computes an actual nonzero Jacobian class D=[P−O]. It checks [373]D and [315]D. Since gcd(373,315)=1, a nonzero class cannot be killed by both. Conditional on the established geometric-quotient/twist dichotomy, these scalar tests distinguish the two signs.

The [certificates](results/quotient-sign-certificates-sage.json) reproduce **eleven positive signs and one negative sign**, with kernel 810 the negative case. Each certificate records its chart and point witnesses, and agrees with the earlier complete Frobenius polynomials. No extension-field curve counts are used to form these certificates. This implements the kind of rationality check described in [Tian, §3.1](https://arxiv.org/html/2012.07173v1#S3.SS1).

This is a conditional test between two known possibilities. It does not classify arbitrary Jacobians from group order alone, compute χ18, or directly evaluate the quotient map and its Frobenius conjugate. The distinction matters because Kummer rationality loses the ± sign. The prior F49 trace-zero witness remains the explicit negative-sign control.

## CM-guided search outside the original example

The search uses three stated cases (p,t,D):

    (13,−83,−211), (13,61,−563), (19,−137,−963).

They satisfy t²−4p³=9D, p³+t+1 divisible by 9, ordinary trace, and prime p³+1−t. Sage computes each Hilbert class polynomial and examines all its roots in the stated cubic field. One model is retained per coefficient-Frobenius j-orbit. An explicit choice between a curve and its quadratic twist gives the requested trace.

These congruences guide the search; the actual module is checked independently. In every retained model, all four roots of the 3-division polynomial lie in the cubic field and all four corresponding y²-values are nonsquares. Frobenius therefore acts as −I on E[3]. Native polynomial-basis arithmetic independently recounts each curve, verifies division roots, and verifies a generator's prime order.

Elements encode a0+p·a1+p²·a2. For p=13, z³+2z+11=0; for p=19, z³+4z+17=0.

| p | A | B | j | Trace | Prime order | Generator |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 13 | 563 | 436 | 47 | −83 | 2281 | (0,257) |
| 13 | 2154 | 1436 | 249 | 61 | 2137 | (1,632) |
| 13 | 28 | 98 | 827 | 61 | 2137 | (0,633) |
| 13 | 1110 | 1428 | 1299 | 61 | 2137 | (3,1) |
| 19 | 1990 | 5893 | 584 | −137 | 6997 | (1,3165) |
| 19 | 3486 | 4648 | 1217 | −137 | 6997 | (3,2891) |

All six have zero rational 2-kernels on the Weil restriction and twelve mixed 3-kernel candidates of the selected type. They give **six distinct Frobenius j-orbits in three elliptic isogeny classes**, with endomorphism fields of discriminants −211,−563,−107. The original class had field discriminant −59. Models with the same p and trace must not be counted as independent elliptic isogeny classes. Their Weil restrictions have the same polynomial as well.

The CM order used to select a root is recorded separately from its endomorphism field. The full endomorphism ring is not claimed to be determined for every new model; −963 is nonfundamental. No new candidate's quotient Jacobian or source-to-Jacobian map is constructed in this batch. They are actual arithmetic candidates, not demonstrated new weak families.

For the first F13³ model, a further bounded experiment follows all four rational 3-isogenies:

| Kernel x | Target A | Target B | Target j | Mixed 3-kernel candidates |
| ---: | ---: | ---: | ---: | ---: |
| 497 | 883 | 1504 | 2139 | 0 |
| 1114 | 1698 | 1486 | 686 | 0 |
| 1229 | 342 | 46 | 682 | 0 |
| 1918 | 286 | 905 | 1906 | 0 |

All four targets retain order 2281, have a single rational 3-division root, and lose the mixed-kernel condition. Explicit Vélu images preserve the generator's order and the scalar-37 relation. This differs from the earlier −59 example's two horizontal and two descending edges: 3 is inert in the field of discriminant −211, so this surface has no horizontal 3-isogenies. The mechanism is established [isogeny-volcano theory](https://arxiv.org/html/1208.5370v3), now checked in another class. These bounded map checks do not enumerate all 2281 group points.

The [Sage report](results/cm-guided-candidates-sage.json) records all class polynomials and selection bounds. The [native report](results/cm-guided-candidates.json) hashes that input and supplies independent arithmetic checks and the four neighbors.

## Wild covers: normalize extension classes, then check the subgroup

The new [Artin–Schreier experiment](artin_schreier_classes.py) handles polynomial covers y²+y=f(x) over odd-degree binary fields. It removes even pole powers by f→f+h²+h, with exact square roots of their coefficients, and normalizes the constant modulo the image of h→h²+h. It returns the actual coordinate correction Y=y+h(x). Multiple poles, rational-function inputs and even-degree constant reduction are outside this implementation.

Sixteen deterministic polynomial presentations over F512 have degree 12 and coefficients outside F8. Each reduces to a genus-2 equation over F8:

    Y²+Y=x⁵+(1+z)x³+x,   z³+z+1=0.

The polynomial identity is checked, and **all 512 affine point transports per presentation** are verified: 8192 maps in total. Reduction commutes with Frobenius in the tested cases. A separate exhaustive oracle over F8 checks all 512 polynomials of degree at most two against all linear-polynomial coboundaries; it obtains sixteen classes of size 32.

This example is a supersingular toy control, and supplies a useful obstruction to overinterpreting the simpler equation. Exhaustive curve counts over F8,F64,F512 are 9,65,513. The complete genus-2 polynomial over F8 is

    P8(T)=T⁴+64,   #Jac(F8)=65.

Cubing its Frobenius roots gives

    P512(T)=T⁴+262144,   #Jac(F512)=262145=65·37·109.

Thus an order-37 or order-109 subgroup in the extension-field Jacobian cannot inject into the rational Jacobian of order 65 over F8. Every homomorphism into that group kills either such subgroup. The change of equation preserves F512-points, and does not change them into F8-points. The remaining factor of order 37·109 belongs in a trace-zero investigation instead. The norm map on divisor classes has not been evaluated here.

This is the next concrete search lesson for wild descent: quotient Artin–Schreier presentations by h²+h before assessing pole orders, then retain the Frobenius factor and subgroup-map checks. The general setting is already studied by [Hess](https://www.cambridge.org/core/journals/lms-journal-of-computation-and-mathematics/article/generalising-the-ghs-attack-on-the-elliptic-curve-discrete-logarithm-problem/315278965D8D277A85812D9498236A88) and [Ballet–Rolland](https://arxiv.org/abs/2505.21656). A genuinely broader search should extend to multiple poles and quotient subextensions; this batch is not such a classification.

## Reproduce and next tasks

From `checker/`, with Sage 10.6 for the three Sage experiments:

```bash
sage -python research/cm_guided_candidates_sage.py --output research/results/cm-guided-candidates-sage.json
python -m research.cm_guided_candidates --output /tmp/cm-guided-candidates.json
sage -python research/quotient_sign_certificates_sage.py --output /tmp/quotient-sign-certificates-sage.json
sage -python research/symmetric_kummer_inversion_sage.py --output /tmp/symmetric-kummer-inversion-sage.json
python -m research.artin_schreier_classes --output /tmp/artin-schreier-classes.json
python -m unittest discover -v
```

The full Python suite has **67 passing tests**, including independent FLINT verification of all recovered divisor supports and residual polynomial identities, using the pinned optional backend. Sage experiments are run separately from Python CI. Reports store input hashes and exact field moduli.

The highest-priority next task is the eight-section calibration and exceptional-divisor handling needed to connect the nonenumerative inverse to the actual quotient map. Next, parameterize reconstruction for the first −211 candidate and compare all four neighbors' routes. Additional targets still requiring work include the opposite-sign trace-zero map, χ18, multiple-pole Artin–Schreier covers, and dimension-four Schottky/Prym tests. None of these pending steps should be inferred from the completed small experiments.
