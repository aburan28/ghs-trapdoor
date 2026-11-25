# Further descent experiments: isogeny routes, twists and target spaces

Two new bounded experiments give concrete leads. An exact 3-isogeny walk connects six curves rejected by the direct mixed-kernel screen to the previously verified descent example. A separate Jacobian calculation finds an order-373 trace-zero divisor on the opposite-sign quartic over F49. Neither calculation establishes a new vulnerable cryptographic family or a practical attack improvement.

The useful search object is a **curve, its reachable isogenous representatives, a polarization, a target and an explicit subgroup map**. Testing one displayed equation or finding a geometric target is insufficient. The priorities below separate completed evidence from experiments still to run.

## 1. Search endomorphism orders and short isogeny routes

This lead directly addresses how one isogenous curve can pass a particular descent test while another fails it. Isogenous curves have the same rational Frobenius polynomial, but their integral torsion modules and endomorphism orders can differ. That can change the available rational polarized quotient kernels. The general mechanism is established isogeny-volcano theory; the new work here is its exact application to the current toy example. See [Sutherland, §§2–3](https://arxiv.org/html/1208.5370v3).

Over K=F7[z]/(z³+z²+1), with encoding a0+7a1+49a2, the earlier source is

    E: y²=x³+248x+42,  #E(K)=373,  j(E)=62.

All four rational cyclic 3-isogenies exist as maps over K, even though none of their nonidentity kernel points are K-rational. Their kernels have x-coordinates 19,65,134,237. The four target j-invariants are 94,320,16,325 respectively.

| Vertices in the complete connected 3-isogeny component | Count | Endomorphism discriminant | Mixed 3-kernel candidates on Res(E) per vertex |
| --- | ---: | ---: | ---: |
| Surface: j=16,62,325 | 3 | −59 | 12 |
| Floor: j=89,94,184,269,291,320 | 6 | −531 | 0 |

All nine models have order 373. The graph has a triangle on the surface, two leaves at each surface vertex, and 18 directed edges in total. Every floor vertex has exactly one outgoing rational 3-isogeny, returning to the surface. The three surface models form one coefficient-Frobenius orbit; the six floor models form two more orbits. These are three j-orbits in one isogeny class, rather than nine independent families. This construction does not certify a practical weakness for any of them.

For an explicit pair, a kernel with x=19 gives

    E' : y²=x³+161x+42,  j(E')=94.

E' has one rational 3-division root, x=13, and its normalized Vélu quotient is the original E. E' has zero mixed kernels of the selected type; E has twelve. Thus **a direct negative result on E' does not exclude the route E'→E→Jac(C)**. The isogeny degrees are coprime to 373, so they preserve that subgroup. This is failure of one construction test, not a curve that intrinsically cannot descend.

The [native experiment](isogeny_descent_routes.py) verifies every rational point and every generator-addition identity on all 18 directed edges: **6,714 point images and 6,714 additions**. It also checks all 373 points on each of the six composed routes into the original source, including the coordinate conjugations. The earlier divisor certificate supplies the final Jacobian map; no new Jacobian inversion is performed. Scalar labels in the report are validation witnesses obtained by tiny-group enumeration, not part of the rational isogeny evaluator. The earlier inversion remains O(r).

The endomorphism explanation is unusually explicit here. With π the K-Frobenius and trace −29,

    disc Z[π] = (−29)²−4·343 = −531 = 3²·(−59).

The surface has π=−1 on E[3], so (π+1)/3 is an integral endomorphism with trace −9 and norm 35. Its discriminant is −59, already fundamental, forcing the maximal order. At a floor vertex this divisibility fails, and the only other possible order has conductor 3 and discriminant −531. Sage independently checks every j against the Hilbert class polynomial for its order; their degrees are 3 and 6.

**Next experiment:** extend the checker to return bounded routes and the least expensive reachable target, with explicit path degree and subgroup preservation. For new cubic-field searches, choose ordinary isogeny classes with 9 dividing t²−4q³, distinguish scalar from Jordan Frobenius on E[3], and walk between orders before doing expensive theta reconstruction. The discriminant congruence alone does not prove scalar Frobenius or a Jacobian quotient.

For the scalar −1 branch, integrality of (π+1)/3 also requires (q³+t+1)/9 to be integral. A more selective parameter generator can enforce this condition before constructing and testing curves. Count distinct Frobenius j-orbits and endomorphism orders, rather than coefficient pairs, and retain the cost of finding a route. No larger-parameter search is completed here.

## 2. Revisit the opposite-sign quartic in a trace-zero subgroup

The kernel-810 quartic from the [reconstruction](QUOTIENT_RECONSTRUCTION.md) has

    P(T)=T⁶−29T³+343,  #Jac(C)(F7)=P(1)=315.

It cannot receive an F7-rational subgroup of order 373. Over F49, however,

    #Jac(C)(F49)=P(1)P(−1)=315·373=117495.

The subgroup ker(Frob7+1) has order P(−1)=373 and lies in Jac(C)(F49). This follows from the degree of Frob7+1 and Frob7²−1=(Frob7−1)(Frob7+1). It is the anti-invariant, or quadratic trace-zero, subgroup. The [Sage experiment](descent_leads_sage.py) constructs an actual nonzero divisor class in it and checks multiplication by 373, not just the order formula.

Write F49=F7[w]/(w²+1). In the original kernel-810 plane model, take

    P=(0:1+3w:1),  Frob7(P)=(0:1+4w:1),
    D=[P−Frob7(P)].

Sage verifies D≠0, 373D=0, and Frob7(D)=−D. Since 373 is prime, D has exactly that order. Of the curve's 50 F49-points, 42 are outside F7 and give nonzero differences of this form. Only the displayed witness's order is separately checked by scalar multiplication. The experiment does not build the quotient map from the elliptic source into this target.

This suggests keeping some geometrically useful targets that fail the original rationality requirement, provided the checker records **the extension and Frobenius eigenspace**. The general twisting issue is established: a geometrically Jacobian abelian threefold can be a quadratic twist of a Jacobian over the base field. See [Ritzenthaler, §§1–2](https://arxiv.org/html/0901.2920v2). Trace-zero groups also have an existing cryptographic literature; see [Lange](https://www.hyperelliptic.org/tanja/preprints/trace_final.pdf) and [Gorla–Massierer](https://www.aimsciences.org/article/doi/10.3934/amc.2015.9.515). Their stated algorithms should not be assumed to transfer unchanged to a nonhyperelliptic genus-3 target.

**Next experiment:** construct the missing equivariant quotient map over F49 and check its image satisfies Frob7(Q)=−Q. Then measure relation generation in that subgroup. A naive factor base using P∈C(Fq²) costs about q² curve points, while the original favorable plane-quartic route uses geometry over Fq. The extension can erase the advantage; existence of the order-373 subgroup alone is not evidence of a speedup.

## 3. Use Frobenius signs to detect misleading Kummer success

There is a precise reason that successful quotient-to-Kummer evaluations do not settle the preceding issue. The Kummer quotient identifies Q with −Q. If Frob(Q)=−Q, its class modulo sign is base-field rational although Q itself is not. Our trace-zero divisor supplies an actual example: its class modulo sign is F7-rational and its nonzero Jacobian class requires F49. This statement follows from the checked Frobenius identity; an explicit coordinate embedding of this witness into Kummer space is not computed here.

**Next experiment:** attach a lift-sign certificate to every Kummer candidate. Compare the quotient map with its coefficient-Frobenius conjugate, then resolve whether the lift obeys +Frobenius or −Frobenius. Use both the matching and opposite-sign quartics as controls. The original kernel-810 obstruction over F7 remains valid.

This can also serve as an early rationality test before enumerating large groups. It must use normalized functions or a chosen differential basis; an arbitrary projective rescaling of theta coordinates can hide the sign being tested.

## 4. Separate the geometry filter from the point-map bottleneck

Two implementation tasks have greater value than another long list of tiny coefficient pairs:

1. **Replace O(r) divisor inversion.** Solve for a degree-3 effective divisor from the six calibrated Jacobian Weil-function values, handling exceptional supports and the remaining global sign. The polynomial-system route is described by [Tian, §3.2](https://arxiv.org/html/2012.07173v1#S3.SS2); [Milio](https://arxiv.org/abs/1709.06063) gives related explicit isogeny methods. Validate on randomly chosen points without scalar labels, then increase field size while keeping the inversion method fixed. Failure to solve the systems reliably is a concrete negative result.
2. **Test a cheaper Serre-obstruction filter.** Compare a correctly normalized χ18 square class with all twelve paired kernel constructions before performing three exhaustive curve counts. The existing 11-to-1 sign split supplies positive and negative controls. Our saved theta values are squared normalized ratios, not automatically a correctly normalized value of χ18: multiplying them without resolving square roots and differential factors would not certify the obstruction. This filter is unimplemented.

Both would improve the research pipeline. Neither, by itself, establishes that a real DLP is cheaper than generic methods. Total comparison must include isogeny search, field extensions, point transfer, relation generation and linear algebra.

## 5. Go beyond cubic elliptic restrictions with polarization checks

The cubic case benefits from the special geometry of dimension three. Moving to a degree-four elliptic restriction, or to a genus-2 Jacobian over a quadratic extension, produces dimension four. An indecomposable principally polarized quotient is then **not automatically a geometric Jacobian**. There is a genuine Jacobian-locus condition; see [Harris–Hulek on the genus-4 Schottky locus](https://arxiv.org/abs/math/0301007).

**Next experiment:** enumerate small rational isotropic kernels on these dimension-four restrictions, retain the polarization data, and test an appropriate Schottky condition before trying to reconstruct a genus-4 curve. Only then check base-field descent and the large subgroup map. This is an experiment design, not an implemented dimension-four checker; reductions of modular-form tests in the chosen finite characteristic need justification.

**Alternative target:** retain a candidate as a Prym rather than forcing it to be a Jacobian. Start with explicit double covers, compute the norm kernel and polarization type, and search for a matching isogeny factor. [Mumford's algebraic Prym theory](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1974a--PrymVar-I-NC.pdf) distinguishes the principal-polarization cases. An arbitrary ramified double cover does not automatically give the principal polarization needed by the existing theta code. The experiment fails as an attack lead if the resulting representation has no cheaper relation algorithm.

## 6. Search wild covers through Artin–Schreier classes

For characteristic-p covers yᵖ−y=f(x), changing f by hᵖ−h gives an equivalent extension after changing y. This suggests searching Frobenius orbits of extension classes, with ramification data, rather than coefficient lists. It is a concrete direction for the user's “wild descent” idea; wild ramification has a technical meaning here and differs from Weil restriction itself.

**Next experiment:** over small composite-degree binary fields, reduce pole terms modulo h²+h, construct Frobenius-stable spans of the resulting classes, and enumerate quotient subextensions. Calculate ramification/genus and full Frobenius polynomials, then test the norm/conorm map on an odd-prime subgroup. A smaller raw coefficient rank without a lower genus or preserved subgroup is not a positive result. Use the earlier checker as a baseline and avoid counting equivalent Artin–Schreier presentations twice.

[Hess's generalization](https://www.cambridge.org/core/journals/lms-journal-of-computation-and-mathematics/article/generalising-the-ghs-attack-on-the-elliptic-curve-discrete-logarithm-problem/315278965D8D277A85812D9498236A88) and [Ballet–Rolland (2025)](https://arxiv.org/abs/2505.21656) are relevant starting points. The latter treats generalized Artin–Schreier descent and explicit small examples, not a new NIST attack. This proposal does not remove the proper-subfield obstruction for prime extension degrees or change the existing NIST binary-curve results.

## Reproduction and evidence boundaries

From `checker/`:

```bash
python -m research.isogeny_descent_routes --output research/results/isogeny-descent-routes.json
sage -python research/descent_leads_sage.py --output research/results/descent-leads-sage.json
python -m unittest tests.test_isogeny_descent_routes -v
```

The native route experiment uses the standard library. The independent Sage 10.6 run checks the 18 Vélu isogenies and all 6,714 point images, the two Hilbert class polynomials, and the trace-zero divisor witness. Its [report](results/descent-leads-sage.json) hashes both input reports. The [native report](results/isogeny-descent-routes.json) also hashes the earlier completed divisor certificate.

Four new Python tests cover complete report reproduction and graph closure, an independent kernel-point-sum oracle over characteristics 5,7,11 with nonrational kernel points, identity and invalid-root cases, and the trace-zero report's scope. The complete Python suite has **59 passing tests**. The Sage computation is separate from Python CI. No results above claim a new cryptographic break, a novel mathematical mechanism, an efficient replacement for divisor inversion, or changed conclusions about NIST curves.
