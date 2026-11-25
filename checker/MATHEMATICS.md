# Mathematical scope and implementation

## Classical invariant

Let `K=GF(2^d)`, `k=GF(2^s)` with `s|d` and `s<d`, `q=2^s`, and `n=d/s`. Let the input ordinary elliptic curve be

`E: y^2+xy=x^3+a*x^2+b`, with `b!=0`.

Form the binary vector space

`U = span_F2 {(1,b), (1,b^q), ..., (1,b^(q^(n-1)))}.`

The classical parameter is `m=dim_F2(U)`. The conventional definition uses square roots of the conjugates of `b`. Squaring is an invertible `F2`-linear map of `K`, commutes with Frobenius, and fixes zero. Consequently the augmented rank and the membership test below are unchanged when using `b` instead of its square root.

The first coordinate is a separate binary coordinate. The rank is **not** taken over `k` or `K`.

Under the classical regularity condition

`n odd OR Tr_(K/F2)(a)=0 OR m=n`,

the genus is

`g=2^(m-1)-epsilon`, where `epsilon=1` if `(1,0)` belongs to `U` and zero otherwise.

The membership criterion avoids extending the odd-degree trace shortcut to inappropriate even-degree cases. For odd `n`, it agrees with `epsilon=1` exactly when `Tr_(K/k)(b)=0`. The checker does not assert the genus result when regularity fails.

The compositum interpretation explains the formula: nonzero characters correspond to vectors in `U`; the characters with both the first and second coordinates nonzero contribute genus one. There are `2^(m-1)` vectors with first coordinate one, with one exception precisely when `(1,0)` is present.

## Subgroup screening

An attack requires the composition of the isogeny transport and the descent map to preserve the target subgroup. A curve/genus calculation alone cannot establish this.

If a claimed subgroup order `r` is provided, the checker applies a necessary Hasse-Weil size condition. It uses the conservative integer upper bound

`#J_C(k) <= (q+1+ceil(2*sqrt(q)))^g`.

For genus one the exact integral upper bound is `q+1+floor(2*sqrt(q))`. Exponentiation is capped at `r` to avoid constructing enormous integers. Failure of this size condition excludes preserving the entire claimed subgroup. Passing the condition gives no information about the actual kernel. No primality, full group order, or nonzero descent image is inferred.

## Explicit degree-3 transport

For this binary model, the 3-division polynomial is

`psi_3(X)=X^4+X^3+b`.

Each root `u` in `K` identifies the Frobenius-stable cyclic kernel `{O,(u,v),(u,v+u)}`. The two nonzero points may be defined over a quadratic extension. Thus limiting the search to `K`-rational points of order three would miss valid rational isogenies.

Vélu's formulas, followed by a y-translation restoring the normalized model, give

`E': y^2+xy=x^3+a*x^2+b'`, with `b'=b+u+u^2`.

For `t=x+u` and points outside the kernel, the map implemented here is

`x' = x + u*x/t^2`,

`y' = y + u*y/t^2 + u*x^2*(x+u+1)/t^3`.

The kernel maps to infinity. Both formulas are evaluated with finite-field operations, and each returned image is checked against the target equation. No auxiliary choice of `v` is needed.

Roots are extracted from `gcd(psi_3,X^(2^d)-X)` and the resulting rational-root factor is split with absolute-trace pairings. Trying polynomial-basis trace functionals separates distinct rational roots deterministically. Small fields use exhaustive roots; tests exercise both methods independently.

Isogenies of degree three are injective on subgroups whose order is coprime to three. This establishes isogeny-path subgroup preservation, **not** preservation by the unconstructed GHS descent map.

This implementation explores only the rational cyclic degree-3 graph. It does not enumerate general prime-degree isogenies, reconstruct Teske's secret path, or exhaust an entire isogeny class. Frobenius/isomorphism invariance is tested, but Frobenius moves are not used as a substitute for new candidate curves.

## Prime-degree field bounds

For an odd prime field degree `d`, write `t=ord_d(2)`. The smallest classical parameter strictly greater than one is `t+1`, yielding the lower genus bound `2^t-1`. These are construction-specific bounds, including all representatives reached by isogenies over the same field when the same construction is used. They do not bound the genus of arbitrary covers or classify other maps out of the Weil restriction.

## Reference fixtures

- **Menezes and Qu**, *Analysis of the Weil Descent Attack of Gaudry, Hess and Smart*, CT-RSA 2001; [preprint](https://cacr.uwaterloo.ca/techreports/2000/corr2000-48.ps), Theorem 1, Section 3, Remark 11, Table 2. This supplies the augmented rank definition, regularity condition, and NIST field-wide bounds.
- **Teske**, *An Elliptic Curve Trapdoor System*, [2003 preprint](https://eprint.iacr.org/2003/058), Sections 1–2 and Appendix A. Both curves use `GF(2)[z]/(z^161+z^18+1)`. With descent to `GF(2^23)`, the supplied secret curve has `m=4, g=8`; the supplied public curve has `m=7, g=64`. The presets retain the polynomial exponents so the transcription can be reviewed.
- **NIST**, [SP 800-186](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186.pdf), Section 3.3, and **FIPS 186-4** for the historical 163-bit specifications. Constants are polynomial-basis parameters. The developer import script cross-checks the constants with OpenSSL's corresponding SEC named parameters. OpenSSL is not a runtime dependency.

Failure of classical GHS is not failure of the general Weil-descent methodology. A useful extension of this checker must provide a concrete alternative construction, prove or verify subgroup preservation, and assess the full search/transport/target-DLP time and memory requirements.
