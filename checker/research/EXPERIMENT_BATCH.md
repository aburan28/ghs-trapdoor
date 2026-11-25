# More experiments: actual torsion modules, Prym factors and tori

This batch adds three reproducible computations to the [alternative-target research map](ALTERNATIVE_TARGETS.md). The parameters are deliberately small enough for exact arithmetic and exhaustive point counts. Positive kernel counts are candidates for quotient construction; no new cryptographic vulnerability is established.

## Cubic-field curve census

The [census](cubic_curve_census.py) checks **16,140 nonsingular short Weierstrass equations**, without deduplicating isomorphic curves. It exhausts all 15,625 coefficient pairs over F125, skipping the 125 singular equations, and samples 512, 64 and 64 nonsingular equations over F343, F1331 and F2197. The samples use the documented SHA256 counter stream without replacement; they do not support population-wide prevalence estimates. Every selected equation is point-counted exhaustively.

| Base field q | Equations checked over Fq³ | Mixed 2-kernel candidates | Mixed 3-kernel candidates | Zero rational 2-kernels, but mixed 3-kernel candidates | Prime-order equations in the preceding column |
| --- | ---: | ---: | ---: | ---: | ---: |
| 5 | 15,500 | 2,542 | 7,688 | 2,604 | 0 |
| 7 | 512 | 81 | 51 | 16 | 2 |
| 11 | 64 | 9 | 39 | 17 | 0 |
| 13 | 64 | 10 | 5 | 1 | 0 |

Each entry counts **equations with at least one kernel candidate**, rather than the total number of kernels. “Mixed” means the necessary no-coordinate-factor condition from the preceding screener. It does not assert geometric indecomposability or a Jacobian over Fq.

The 2-division cubic has 0, 1 or 3 rational roots, determining the E[2] Frobenius conjugacy class. The 3-torsion calculation uses

\[
\psi_3(x)=3x^4+6Ax^2+12Bx-A^2.
\]

For distinct Frobenius eigenvalues, trace and determinant determine the conjugacy class. For a repeated eigenvalue, a scalar matrix fixes all four projective torsion lines, whereas a nontrivial Jordan block fixes one. Counting the rational roots of ψ₃ therefore resolves this ambiguity. The output gives matrices **up to conjugacy**, not a coordinate basis of E[3]. Tests independently enumerate the rational order-3 points on the F343 representatives and their quadratic twists using the group law, checking both eigenspaces.

The two prime-order candidates are over K=F7[z]/(z³+z²+1), with integer encoding a₀+7a₁+49a₂ for a₀+a₁z+a₂z²:

| Equation y²=x³+A·x+B | #E(K) | j | A verified order-373 point |
| --- | ---: | ---: | --- |
| A=248=3+5z², B=42=6z | 373 | 62=6+z+z² | (0,176) |
| A=3, B=76=6+3z+z² | 373 | 16=2+2z | (0,119) |

Neither j lies in F7. Each has no rational (2,2,2) kernel on its cubic Weil restriction, while the 3-torsion module has 16 rational maximal-isotropic kernels, 12 passing the no-coordinate-factor test. Frobenius on E[3] is −I: all four ψ₃ roots lie in K and there are no nonzero rational order-3 points. The group-order and point-order checks are exact at these toy sizes.

This is additional input data for the mechanism in [Tian's cover-attack paper](https://arxiv.org/html/2012.07173v1), which already treats prime-order cubic-extension curves. The computation does not establish that these particular kernels yield a suitable Jacobian over F7. Target reconstruction, its rationality obstruction, explicit subgroup transport and an attack-cost comparison remain open steps.

A useful search filter follows from the module survey: at ell=3 and q≡2 mod 3, a mixed candidate requires trace≡0, hence #E(Fq³)≡0 mod 3. This excludes prime group orders greater than 3 for that case. At q≡1 mod 3 the prime-order candidates in this test require scalar −I; scalar +I forces rational 3-torsion. This is a filter for the stated kernel construction, not for all Weil descent.

The [saved report](results/cubic-curve-census.json) includes field moduli, aggregate module classes, exact representative equations, all prime-order extension-j candidates encountered, subgroup witnesses, bounds and SHA256 digests of the ordered per-equation data. It does not store all 16,140 equation records; rerunning reconstructs the stream and its digests.

## An explicit Prym target

For a squarefree monic quartic f with f(0)≠0, the [cover experiment](prym_cover_experiments.py) constructs

\[
C:y^2=f(x^2),\qquad E:y^2=f(u),\qquad D:v^2=u f(u),
\]

of genera 3, 1 and 2. The two degree-2 maps are

\[
\pi_E(x,y)=(x^2,y),\qquad \pi_D(x,y)=(x^2,xy).
\]

These are standard involution/quotient constructions; see the discussion of genus-3 involutions and complementary quotients in [Borówka–Shatsila, §3](https://arxiv.org/html/2308.07038#S3). That paper studies a more specialized family over the complex numbers. The finite-field claims here follow from the displayed maps, the usual norm/pullback identities and the computed Frobenius polynomials, rather than an assertion that its moduli results automatically apply in these characteristics.

On the pullback image of Jac(D), the involution (x,y)↦(−x,y) acts as −1. The norm to E therefore vanishes there, placing that image in Prym(C/E), of dimension two. The identity πD,*πD*=[2] shows that every odd-order subgroup of Jac(D) survives pullback. This is a theorem about the maps; pointwise Jacobian divisor arithmetic is not implemented. The genus-2 quotient map C→D is unramified and can have a 2-torsion pullback kernel, so injectivity is not asserted on all torsion.

The batch takes the first eight admissible quartics in the specified coefficient order for each p=5,7,11,13. All **32** cases verify the full characteristic-polynomial identity

\[
P_C(T)=P_E(T)P_D(T)
\]

from exact point counts over Fp, Fp² and Fp³. Both points at infinity on the monic even-degree models are included. In **23/32** cases, the order of Jac(D)(Fp) contains an odd prime absent from #E(Fp); its primary component belongs to the Prym factor rather than the elliptic factor.

For example, at p=7 take f(u)=u⁴+6u³+1. Then

\[
P_E(T)=T^2+4T+7,\quad
P_D(T)=T^4+4T^3+10T^2+28T+49,
\]

and #E(F7)=12, #Jac(D)(F7)=92=4·23, #Jac(C)(F7)=1104. The order-23 component survives pullback and is absent from the elliptic factor. This gives an explicit experiment with a Prym target and verified isogeny factors. It does not provide a reduction from an arbitrary elliptic DLP or demonstrate a faster DLP algorithm.

Equations, all point counts, full polynomials and map-check counts are in the [Prym report](results/prym-cover-experiments.json). Affine maps are checked pointwise; exceptional infinity charts are not separately evaluated.

## Torus transfer and smooth negative controls

The [torus experiment](torus_transfer_experiments.py) uses the established singular-curve model described by [Kohel](https://www.i2m.univ-amu.fr/perso/david.kohel/pub/torus.pdf). For nonsquare δ∈Fp, the smooth locus of the nodal curve

\[
y^2=x^2(x+\delta)
\]

maps to the norm-one torus in Fp² by t=(y−γx)/(y+γx), γ²=δ. The identity maps to 1; the node (0,0) is excluded. The inverse uses z=γ(1+t)/(1−t), x=z²−δ, y=zx.

At p=5,7,11,13, the script verifies the entire bijection, its inverse and every pair in the group law: **36+64+144+196=440 pairs**. It transfers and solves all **40** possible scalar logs on these toy cyclic groups by exhaustive enumeration in the torus. This is a correctness control, not a cryptographic-size speed comparison.

For each p it then tests the same candidate rational function on the smooth elliptic perturbation y²=x³+δx²+1. All four controls produce an explicit pair violating the homomorphism law. Thus the singular-model construction does not automatically extend to that smooth curve. This refutes only this candidate formula; it is not a test of every possible reduction. See the [torus report](results/torus-transfer-experiments.json).

## Reproduce and extend

Run from `checker/` with Python 3.10 or newer; no external algebra system is needed:

```bash
python -m research.cubic_curve_census --output /tmp/cubic-curve-census.json
python -m research.prym_cover_experiments --output /tmp/prym-cover-experiments.json
python -m research.torus_transfer_experiments --output /tmp/torus-transfer-experiments.json
python -m unittest discover -v
```

The follow-up [quotient reconstruction](QUOTIENT_RECONSTRUCTION.md) now constructs all 24 target curves, records two opposite-sign obstructions, and verifies actual subgroup maps for both Galois-related input models. Its divisor inversion enumerates the tiny target group. The next priority is efficient polynomial-system inversion and a complete cost comparison against Pollard rho and known cubic-field index calculus. Other unrun experiments are quartic-extension torsion kernels, Prym divisor transport, and longer genus-2 isogeny paths. Those remain proposals, rather than results of this batch.
