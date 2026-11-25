# Explicit quotients for the two cubic-field candidates

The two prime-order inputs from the [curve census](EXPERIMENT_BATCH.md) now have explicit genus-3 target curves. This is a toy reconstruction of the established theta-quotient method, with exact finite-field arithmetic and a separate experiment for Jacobian divisor transport. It does not establish a new cryptographic vulnerability or a practical speedup.

## Curves and the rationality distinction

Both source curves are over K=F7[z]/(z³+z²+1), use integer encoding a₀+7a₁+49a₂, and have #E(K)=373:

| Source y²=x³+A·x+B | A | B | Verified order-373 generator |
| --- | ---: | ---: | --- |
| First census candidate | 248 | 42 | (0,176) |
| Second census candidate | 3 | 76 | (0,119) |

The two inputs are **Galois-related coordinate instances of one Weil-restriction case**. Their j-invariants have Frobenius orbit [62,16,325]. With u=78=1+4z+z², the map (x,y)↦(u²x⁷,u³y⁷) transports the first subgroup to the second and sends the displayed generators to each other; all 373 points are checked. The second curve is K-isomorphic to the coefficient-conjugate of the first, and their Weil restrictions over F7 are isomorphic.

An explicit 6×6 torsion-module isomorphism commutes with Frobenius and pairs all twelve kernels across the two source models. Thus the 24 constructions give **at most twelve** distinct quotient classes. Further deduplication within those twelve kernels has not been performed. The native Kummer reports select kernel 153 on each source model. The completed Sage divisor transfer uses kernel 153 of the first source and its corresponding kernel 397 on the second, sharing one verified Jacobian target.

The earlier screener gave 12 rational maximal-isotropic 3-torsion kernels per source after the no-coordinate-factor test. The [new implementation](quotient_reconstruction.py) reconstructs a smooth plane quartic over F7 for **all 24** kernel constructions, and verifies all 28 bitangent lines on each curve. These counts are of kernel constructions; targets have not been deduplicated by isomorphism class.

| Each source's 12 kernels | Count | Curve counts over F7,F49,F343 | Frobenius polynomial | Jacobian order over F7 |
| --- | ---: | --- | --- | ---: |
| Matching target class | 11 | 8,50,431 | T⁶+29T³+343 | 373 |
| Opposite-sign target class | 1 | 8,50,257 | T⁶−29T³+343 | 315 |

The Weil restriction has characteristic polynomial T⁶+29T³+343. The computations use three exhaustive projective point counts, which determine the **complete** genus-3 polynomial, rather than just comparing group orders. All 24 quartics have exact smoothness certificates from affine resultants, explicit treatment of extraneous resultant roots, and an infinity-chart check.

The opposite-sign construction is a concrete rationality obstruction: a group of order 315 cannot contain the source's order-373 subgroup. Its polynomial is the quadratic-twist polynomial of the required abelian variety. A geometric Jacobian reconstruction therefore does not automatically give the required Jacobian over F7. This is the familiar genus-3 twisting/Serre-obstruction issue, not a failure of Weil restriction or proof that every kernel fails. The other eleven constructions per source have the required full polynomial.

For the first source and kernel index 153, one matching quartic is

\[
\begin{aligned}
0={}&5X^4+X^3Y+5X^3Z+6X^2YZ+6X^2Z^2\\
 &+6XY^3+6XY^2Z+XYZ^2+4XZ^3\\
 &+2Y^4+6Y^3Z+6Y^2Z^2+YZ^3.
\end{aligned}
\]

The [saved report](results/quotient-reconstruction.json) contains all kernel bases, seeds, quartic coefficients, counts, polynomials and map witnesses. Kernel indices refer to the exhaustive Lagrangian enumeration; the actual basis is stored so the index is not needed to interpret the subgroup.

## Construction and checks

The implementation follows the normal-Weil-function, theta and Riemann-quartic approach in [Tian, §§2–3](https://arxiv.org/html/2012.07173v1). Tangent-line functions give the elliptic 3-Weil functions, whose normalization identities are checked on all 81 torsion pairs. Product functions and the rational kernel then define an invariant theta section. This is additional explicit data for a published mechanism.

The arithmetic backend is python-flint 0.9.0. A fixed, irreducibility-checked degree-18 modulus supplies the 2-division splitting field and the 3-torsion points; every extension and square root is exact. Root and square-root choices are made canonical. The F343 embedding is independently checked against the existing polynomial-basis implementation.

There is a small but consequential implementation trap: projecting the constant seed onto the invariant section can give the **zero section**. This occurs in four of the 24 constructions. It is not a decomposability certificate. The implementation tries the even tensor seeds in the basis {1,x,y}⊗³ of L(3Θ) and recovers a nonzero section in each case. All recovered sections give 36 nonzero even theta-null values and 28 nonzero odd-theta gradients.

Those first jets give the 28 bitangents in the source Weil restriction's tangent coordinates. Reconstructing from an Aronhold system in a rational basis of Lie(Res E) yields equations directly over F7. Two independent Aronhold systems reproduce the same normalized quartic in a regression check. The use of gradients to obtain bitangents is classical; see [Fiorentino, §4](https://arxiv.org/abs/1612.02049). The Lie-basis calculation is semilinear Galois descent, and is distinct from differentiating the inseparable Frobenius endomorphism.

The 64 homogeneous quotient-to-Kummer coordinates are evaluations of

\[
g_T(R)^3\,\theta(R+T)^2,\qquad T\in W[2].
\]

They avoid division by θ(R), including its zero locus. The code evaluates these sections on all 373 points of each source subgroup for one matching kernel per source. Each gives exactly **187** projective images, with only P and −P identified. At each of three source points, translating by all 27 kernel elements preserves the image: **81** checks per source. A translation outside the kernel changes the image in a negative control. Exceptional poles outside these evaluated points and translates are not covered by this bounded evaluator.

## Jacobian divisor experiment

The separate [Sage script](quotient_divisor_verification.py) uses actual genus-3 Jacobian arithmetic, not an arbitrary isomorphism between equal-order finite groups. The [exported input](results/quotient-divisor-inputs.json) gives source points, 64 theta-null values, the 28 labeled bitangents, and the seven selected quotient-section values for every source point. Scalar labels are supplied for validation; they are not a source-DLP oracle used to evaluate the theta functions.

The bitangent contact divisors supply 2-torsion divisors and the even theta characteristic. Riemann-Roch bases construct six Jacobian Weil functions. Common local uniformizer factors cancel the basis poles at support points, and the corresponding factors cancel in the norm terms. The resulting Weil functions are calibrated at the identity to the source theta-null values. This links the explicit quotient coordinates to divisor classes on the reconstructed curve. Sage's [Hess model](https://doc.sagemath.org/html/en/reference/function_fields/sage/rings/function_field/jacobian_hess.html) supplies the target group arithmetic.

Divisor inversion in this experiment **enumerates the tiny order-373 Jacobian group modulo sign**. An additional quotient evaluation at P+G resolves the sign using the chosen image of G. Its O(r) inversion cost is infeasible at cryptographic subgroup sizes. Even a successful full toy transfer is therefore a correctness result, not an attack-cost improvement. The polynomial-system inversion described in the published method remains a separate implementation task.

The [completed Sage report](results/quotient-divisor-verification.json) records **373 lifted points, 373 generator-addition checks, and the scalar-37 relation for each source model**. The first lift evaluates actual Jacobian Weil functions on all 187 classes modulo sign. The second reuses that certified lift through the independently checked conjugation/scaling map, including its kernel correspondence. The addition checks use the rational Jacobian group over F7, with an additional scalar check in the extension-field group used for the Weil functions.

For a concrete divisor witness, write F for the first quartic displayed above and use C':F(X,Y+Z,Z)=0. Let O=(1:4:1), P0=(1:1:0), and Q=(1:2:4). With a consistent choice of the global sign, both displayed source generators map to

\[
110[P_0-O]=[Q+2P_0-3O].
\]

The report also gives the degree-3 closed-place divisor for the scalar-37 image. Thus this checks actual subgroup maps into a stated Jacobian, beyond equality of Frobenius polynomials or an arbitrary equal-order group identification. The inversion remains O(r). The native Python reconstruction report keeps its divisor-lift flag false because that verification is performed by the separate Sage experiment.

The default Sage command handles the Galois-related pair through this shared target. The optional `--source-index 1` command instead requests a separate inversion for kernel 153 of the second model; that additional target is outside the completed default report.

## Reproduce

From `checker/`:

```bash
python -m pip install -r research/requirements-reconstruction.txt
python -m research.quotient_reconstruction --output /tmp/quotient-reconstruction.json
python -m research.quotient_reconstruction --export-divisor-inputs --output /tmp/quotient-divisor-inputs.json
python -m unittest discover -v
sage -python research/quotient_divisor_verification.py --input /tmp/quotient-divisor-inputs.json --output /tmp/quotient-divisor-verification.json
```

The complete Python suite has **55 passing tests** with the optional backend installed. Ten reconstruction tests cover independent field arithmetic and cocycle probes, alternate Aronhold systems, zero-seed recovery, source-equivalence and kernel-pairing checks, all saved reconstruction and subgroup data, smoothness controls, kernel/nonkernel translations and rejection of optimized Python that would strip assertions. The CI workflow installs the pinned optional backend on Python 3.10 and 3.13; the binary checker itself keeps its dependency-free runtime.

Sage 10.6 is used for the separate divisor experiment. It is not part of the Python CI job. Native source installation, or the official `sagemath/sagemath:10.6` container, can run that script. The tested container digest is `sha256:19995db6194f4a4bab18ce9a88556fd15b9ed5e916b4504fefe618a7796ddbdb`.
