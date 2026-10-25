# Genus-2 isogeny and Weil-descent audit

The [calibration experiment design](CALIBRATION_EXPERIMENT.md) specifies the next test, with a completed source-side rank/frame preflight and 32 frozen challenge inputs. It separates calibration using the previous toy anchor from initialization without that anchor. Target calibration and blind transfer have not run.

The latest [follow-up experiments](DESCENT_FOLLOWUP.md) add three nonenumerative algebraic inverse fibers in a symmetric Kummer model, conditional sign certificates, six CM-guided candidates in three additional isogeny classes, and sixteen Artin–Schreier presentation controls. The complete suite has **67 passing Python tests**; three new Sage experiments also completed. Connecting the algebraic inverse to the actual elliptic quotient still requires coordinate calibration.

The [further descent leads](DESCENT_LEADS.md) add an exact nine-vertex 3-isogeny graph: six curves with no direct mixed 3-kernel route back to the earlier example in one step. A separate Sage computation constructs an order-373 trace-zero Jacobian divisor on the opposite-sign target over F49. The notes give six research directions with concrete next experiments and evidence limits. At the #7 stage, the suite had **59 passing Python tests**, plus the completed independent Sage experiment.

The [new experiment batch](EXPERIMENT_BATCH.md) adds an actual cubic-field curve census, explicit Prym covers and torus transfer controls. The [quotient reconstruction](QUOTIENT_RECONSTRUCTION.md) now constructs 24 plane quartics from the two prime-order F343 inputs, distinguishes 22 matching targets from two opposite-sign cases, and evaluates explicit quotient-to-Kummer maps. The #6-stage source-tree suite had 55 passing tests with the optional reconstruction backend installed. Reports include exact parameters and reproducible computations.

The published GLV/GLS lead produced a negative result for one precise route, plus a reproducible toy example of why isogenies can change a direct descent test. **No previously unreported cryptographic vulnerability has been established.**

The inputs are the exact parameters in Appendices A.2 and A.3 of Bos–Costello–Hisil–Lauter's [March 2013 preprint](https://www.microsoft.com/en-us/research/wp-content/uploads/2013/03/8glv.pdf). The generic and Kummer entries use the same three underlying curves. The authors already considered Weil descent and index calculus; this audit does not treat them as previously believed unaffected by those methods. “NIST-friendly” below describes a prime chosen for arithmetic, not a NIST standardized curve.

| Published source | Base prime p | Rational Richelot kernels | Immediate quotients passing the branch test |
| --- | --- | --- | --- |
| BK, Montgomery | (2^31−201)·2^32−1 | 0 | No first step |
| BK, NIST-friendly | 2^64−2285 | 0 | No first step |
| BK, Mersenne | 2^61−1 | 0 | No first step |
| Generic/Kummer, Montgomery | (2^31−307656)·2^32−1 | 15 | 0/15 |
| Generic/Kummer, NIST-friendly | 2^64−189 | 15 | 0/15 |
| Generic/Kummer, Mersenne | 2^61−1 | 15 | 0/15 |

All three generic/Kummer source curves also fail the branch test. All 45 quotient constructions have nonzero Richelot determinant. These are counts of kernel/quotient constructions, not a claim that 45 different isomorphism classes were obtained. The saved [report](results/genus2-audit.json) contains the parameters, all partitions, determinants, quadratic factors and codomain equations.

The BK obstruction is exact for this isogeny type over K=Fp². Since fifth roots of unity lie in K, the computed nontrivial fifth-power residue of −a proves x^5+a irreducible. K-Frobenius therefore fixes infinity and cycles the other five branch points. A rational Richelot kernel would require an invariant partition into three pairs; the partner of infinity would have to be another rational branch point, which does not exist. This rules out a rational Richelot first step, not arbitrary odd-degree isogenies or constructions over larger fields.

The selected target class is

\[
C/K:\quad y^2=\epsilon(x-\alpha)h(x),\qquad
\alpha,\epsilon\in K,\quad h\in\mathbb F_p[x],\quad \deg h=4\text{ or }5.
\]

This is the quadratic-extension genus-2 class in §4.2.1 of [Momose–Chao](https://eprint.iacr.org/2005/277). We test a necessary condition under every change of x-coordinate in PGL₂(K), rather than testing the displayed coefficients alone. A negative result excludes that class; it does not exclude other covers or generalized Weil descent. A positive result is recorded as branch data requiring further checks, not as an attack verdict.

Write σ(x)=x^p on the geometric projective line. In any model of the selected form, one K-rational exceptional branch can be removed, leaving five branches invariant under a conjugate τ=Mσ of Frobenius. Necessarily M∈PGL₂(K) and Mσ(M)=1. For each possible exceptional branch, the code enumerates all 5·4·3 possible images of three of the remaining branches. Three point images determine M uniquely, so at most 360 possibilities cover all relevant changes of coordinates. Roots of the quotient quadratics are represented in Fp⁴; the code checks that the resulting M is actually defined over K. This is an exhaustive test of this branch condition, not of all possible descents.

Run the published audit from the project root:

```bash
python research/genus2_audit.py --output research/results/genus2-audit.json
```

Custom Rosenhain inputs describe y²=x(x−1)(x−λ)(x−μ)(x−ν). Each coefficient is `[real, imaginary]` in Fp[i]/(i²+1). The implementation accepts primes p≡3 mod 4 below 2^64, and checks primality and nonsingularity. The source Rosenhain model has all six branch points rational over K, so all 15 first-step kernels are rational.

```json
{"family":"rosenhain","p":7,"parameters":[[2,0],[3,0],[1,1]]}
```

```bash
python research/genus2_audit.py --input research/examples/weak-branch-toy.json
```

A BK input instead uses `{"family":"bk","p":11,"a":[1,1]}`. This test additionally requires 5 to divide p²−1. The BK command reports the first-step kernel obstruction; it does not enumerate quotient models for a split BK source. The default published inputs all have irreducible quintics.

The explicit toy construction uses K=F7[i]/(i²+1):

\[
C:\quad y^2=(x-(1+i))x(x-1)(x-2)(x-3).
\]

Partition the branches `[infinity,0,1,2,3,1+i]` as `((0,1),(2,4),(3,5))`. The resulting Richelot quotient is

\[
D:\quad Y^2=(5+4i)X^6+2X^5+6X^4+(4+3i)X^3
             +5iX^2+(5+5i)X+3i.
\]

D fails the selected branch condition under every PGL₂(K) change of coordinate. Its dual Richelot isogeny maps Jac(D) back to Jac(C), where the known cover construction applies. Thus failure of the direct normal-form test need not survive an isogeny.

An actual smooth genus-4 cover over F7 is

\[
H:\quad S^2=2t^9+t^8+5t^7+5t^3+6t^2+2t,
\]

with map to C over K

\[
x=\frac{-t^2+2t+1}{2t},\qquad
y=\frac{(1-it)S}{8t^3}.
\]

The script verifies squarefreeness and the rational-function identity coefficient by coefficient, and checks all 40 affine H(F49) points with t≠0. It does not separately evaluate the exceptional coordinate charts at zero and infinity. Exhaustive point counts give C and D counts 46 and 2334 over F49 and F2401, while H has 8,42,344,2266 points over F7,F49,F343,F2401. Consequently their complete Frobenius polynomials satisfy

\[
P_C(T)=P_D(T)=T^4-4T^3-26T^2-196T+2401,
\qquad P_H(T)=P_C(T^2).
\]

This verifies the Weil-restriction isogeny class, using more than equality of group orders. All three relevant rational Jacobian groups have order 2176=2^7·17. The degree-4 Richelot isogeny and its dual preserve the order-17 subgroup. Pointwise transport of Jacobian divisors and a DLP computation are not implemented here. These parameters are too small to claim a meaningful attack speedup; the example illustrates a known mechanism, not a new vulnerable family.

```bash
python -m research.genus2_cover_example --output research/results/toy-cover.json
python -m unittest discover -v
```

The genus-2 audit adds nine experimental checks to the original 16 binary-checker tests. The [alternative-target extension](ALTERNATIVE_TARGETS.md) adds nine more checks, the [new batch](EXPERIMENT_BATCH.md) adds eleven, and the [quotient reconstruction](QUOTIENT_RECONSTRUCTION.md) adds ten, bringing the full source-tree suite to 55 tests. Independent validation includes an exhaustive PGL₂(F9) oracle on every set of six rational branches, an oracle with nonsplit branches, coordinate changes on a known positive example, and full Frobenius-polynomial comparisons for all 15 Richelot quotients of a separate toy curve. The source tree includes this research directory; installing the binary-checker wheel does not install it.

The remaining research leads are hypotheses, with concrete conditions that the present audit does not establish:

- Longer paths from the three generic/Kummer sources could reach a different weak representative. Any extension of the search must enumerate Frobenius-stable kernels even when individual kernel points or quadratic factors are not rational; testing only split factorizations would miss possibilities.
- Isogenies of the Weil restriction itself can reach targets inaccessible to an elliptic-curve walk. [Tian's cubic-extension construction](https://arxiv.org/html/2012.07173v1) is a published example of this approach, not a new finding. For the quadratic genus-2 case the Weil restriction has dimension four, and a candidate quotient must actually be a suitable Jacobian, with its rationality, polarization, subgroup map and cost checked.
- Alternative covers beyond the selected normal form remain untested. The [2009/2013 classification](https://eprint.iacr.org/2009/613) concerns its stated (2,…,2)-cover setting and does not classify all possible reductions. The [2025 Artin–Schreier descent paper](https://arxiv.org/abs/2505.21656) supplies further descent theory and small examples; it does not establish a new practical attack against the audited curves.

No conclusion here changes the NIST binary results in the main checker. In particular, neither a negative normal-form test nor the absence of rational Richelot kernels is a universal security certificate.
