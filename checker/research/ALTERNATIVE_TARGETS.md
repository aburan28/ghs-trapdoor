# Weil descent beyond hyperelliptic curves

The research question is broader than finding another GHS cover: can the same scalar relation be transferred efficiently to a group, or a compatible group action, with cheaper discrete-log algorithms? The term is **Weil descent**. There is no general theorem that a failed classical GHS test rules out every such transfer.

For a curve or Jacobian A over K=Fq^n, Weil restriction gives W=Res(K/Fq)(A), with W(Fq) canonically equal to A(K). It changes the field and dimension, not the group size. A useful attack still needs a computable transfer preserving the relevant subgroup and a net improvement after map construction, relation generation, linear algebra and individual logarithms. For an isogeny φ of degree coprime to the subgroup order r, injectivity on that subgroup is automatic. For other maps it must be proved separately.

The literature already provides several different spaces. None of the connections in this table is claimed as a new discovery.

| Space or curve type | Established connection | What an additional attack would need |
| --- | --- | --- |
| Weil restriction itself, as an abelian variety | Gaudry's index calculus works directly on suitable abelian varieties, without constructing a GHS covering curve. [Paper](https://www.sciencedirect.com/science/article/pii/S074771710800182X). | A demonstrable improvement over that existing baseline, with realistic relation-solving costs. |
| Non-hyperelliptic Jacobians | Tian constructs suitable quotients of cubic elliptic Weil restrictions as genus-3 Jacobians, including a prime-order example. [Paper](https://arxiv.org/html/2012.07173v1). | A rational polarized kernel, a suitable Jacobian over Fq, an explicit point map and competitive cost. The new experiment below tests the kernel stage. |
| Prym varieties | For a curve cover, the connected kernel of the Jacobian norm is an abelian variety, complementary up to isogeny to the lower curve's Jacobian. [Finite-field Prym paper](https://arxiv.org/abs/0706.0121). | A cover/correspondence placing the intended r-subgroup in that factor, plus an efficient DLP algorithm on its representation. This is an untested attack lead here, not a reported break. |
| Trace-zero varieties | Kernels of trace maps admit their own index-calculus analysis. [Gorla–Massierer](https://arxiv.org/html/1405.1059v2). | The correct source subgroup: the standard construction starts from a curve defined over the base field. It does not automatically compress a generic E/Fq^n into dimension n−1. |
| Tori and generalized Jacobians of singular curves | Kohel realizes algebraic tori through singular-curve divisor arithmetic and transfers certain finite-field multiplicative DLPs. [Paper](https://www.i2m.univ-amu.fr/perso/david.kohel/pub/torus.pdf). Generalized-Jacobian DLPs also decompose through their affine and abelian factors: [Galbraith–Smith](https://arxiv.org/html/math/0610073). | An explicit faithful map from the original source. Introducing singularities changes the group; changing the equation to a singular curve does not supply such a map. |
| Kummer surfaces and other geometric quotients | Generalized Kummer surfaces can give compact models for quadratic-extension elliptic points. [Koshelev](https://doi.org/10.1016/j.ffa.2020.101774). | A cheaper way to recover the scalar. A rational parameterization or efficient point compression alone does not linearize the inherited group action. |
| Superelliptic and higher-genus Jacobians | Generalized descent already extends beyond the original binary elliptic setting. [Iijima–Shimura–Chao–Tsujii](https://doi.org/10.1093/ietfec/e88-a.1.97), [Hess](https://www.cambridge.org/core/journals/lms-journal-of-computation-and-mathematics/article/generalising-the-ghs-attack-on-the-elliptic-curve-discrete-logarithm-problem/315278965D8D277A85812D9498236A88). | Family-specific cover genus, subgroup preservation and cost analysis. Merely choosing a new curve equation is not evidence of a new reduction. |
| p-adic/formal-group lifts | Finite residue fields connect to unramified p-adic extensions and Witt vectors. Smart's trace-one attack is a special useful lift: [original report abstract](https://shiftleft.com/mirrors/www.hpl.hp.com/techreports/97/HPL-97-128.html). | A computable relation-preserving lift with nonzero useful logarithms and controlled precision. Arbitrary coordinate lifts do not preserve Q=[x]P. |

There are also strong, established links from finite-field geometry to number fields: Honda–Tate organizes simple isogeny classes by conjugacy classes of Weil numbers. [Honda's original paper](https://www.jstage.jst.go.jp/article/jmath1948/20/1-2/20_1-2_83/_article) and [Milne's early account](https://www.jmilne.org/math/articles/1971a.pdf) are relevant starting points. This makes Frobenius and endomorphism data useful for organizing target searches. It does not itself compute pointwise isogenies or discrete logarithms.

A field embedding from a characteristic-p field into a characteristic-zero field cannot preserve 1 and addition: applying it to p·1=0 would give p·1=0 in characteristic zero. Witt/Teichmüller lifts are therefore richer lifting constructions, not ordinary field embeddings. For multiplicative finite-field points the Teichmüller lift preserves multiplication, but the p-adic logarithm kills these roots of unity. More generally, a homomorphism from an r-torsion group to a characteristic-zero additive group sends every element to zero, since r·log(P)=log([r]P)=0. [Sage's documented logarithm implementation](https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/padic_generic_element.html) states and tests this torsion obstruction. The anomalous-curve attack avoids a naive torsion embedding by its special lifting construction; it is not a general consequence of moving to Qp.

Likewise, properness implies that a global algebraic group morphism from an abelian variety to an affine algebraic group is constant. Pairing reductions use a different construction; one cannot substitute arbitrary coordinate functions for a pairing and assume they give a faithful group homomorphism. An abstract embedding of a finite cyclic group into roots of unity is also insufficient unless it can be evaluated from the input points without first knowing the discrete log.

## Experiment: rational quotients of a cubic Weil restriction

The new [kernel screener](weil_quotient_kernels.py) works on W[ell], for W=Res(Fq^3/Fq)(E), at ell=2 or 3 coprime to the characteristic. It requires the actual 2×2 action M of q^3-Frobenius on E[ell]. In bases on the conjugate factors obtained by q-Frobenius, the alternating pairing has block weights 1,q,q². The 6×6 q-Frobenius cycles the factors and applies M at the wraparound. Its pairing multiplier is q, and its cube acts by M on every factor.

The code enumerates every dimension-3 isotropic subspace L of F_ell^6, tests Frobenius stability, then checks L∩E_i[ell]=0 for each geometric elliptic factor. There are exactly (ell+1)(ell²+1)(ell³+1) geometric Lagrangians: 135 for ell=2 and 1120 for ell=3. This exhausts the relevant prime-level kernels for the supplied torsion module. The no-coordinate-factor test is necessary for the intended gluing route, not a proof of geometric indecomposability or of a Jacobian model over Fq.

For the actual example in Tian §3.4, K=F523[w]/(w³+2),

\[
E:\ y^2=x^3+(82w^2+92w+140)x+359w^2+339w+243.
\]

We verify two independent 3-torsion points from the paper: curve membership, the 3-division polynomial, and distinct x-coordinates. Their x-coordinates lie in K and their y-coordinates are multiples of θ, where θ²=w+6 is a nonsquare in K. Thus q^3-Frobenius is −I on the entire E[3], rather than a matrix guessed from the point count. Independently, a polynomial gcd with X^(523³)−X proves the 2-division cubic irreducible, so 2-Frobenius has the irreducible order-3 action.

| Same published curve | Geometric kernels | Rational kernels | Passing the no-coordinate-factor test |
| --- | --- | --- | --- |
| ell=2 | 135 | 0 | 0 |
| ell=3 | 1120 | 16 | 12 |

The twelve entries are kernel subgroups, not a claim of twelve nonisomorphic curves or attacks. The supplied subgroup order 143069629 is from the paper; our code does not recount E(K). It is coprime to 27, so a quotient isogeny of that degree would preserve a subgroup of that order if the claimed input is valid. The paper already constructs one target plane quartic and transfer. We do not reconstruct its theta functions, quotient equations or point map here. This is a reproduction and extension of the screening stage of a known attack, not a new vulnerable family.

A useful failure mode appears in the finite module survey: M=−I and M=[[0,1],[2,1]] over F3 have the same trace and determinant, but give 12 versus zero kernels passing the mixing test. The nonscalar case has just one rational Lagrangian. Therefore a screener using only the group order or characteristic polynomial can give the wrong kernel answer when Frobenius has repeated eigenvalues. With distinct roots the characteristic polynomial fixes the 2×2 conjugacy class; the repeated-root case needs more torsion data.

Run from the repository root:

```bash
cd checker
python -m research.weil_quotient_kernels --output /tmp/weil-quotient-kernels.json
python -m unittest discover -v
```

Custom input is `{"q":523,"ell":3,"frobenius_matrix":[[2,0],[0,2]],"subgroup_order":143069629}` passed with `--input FILE`. The 2×2 matrix uses row-vector convention. Custom q is assumed to be a prime power, and the input matrix's realization on the claimed curve is not verified. Built-in curve realization is checked using the torsion calculations above. Output records these distinctions.

The [saved report](results/weil-quotient-kernels.json) includes explicit kernel bases, induced Frobenius actions, subgroup conditions and a complete residue-class module survey. The survey is a finite linear-algebra classification, not a count of vulnerable curves. Nine added tests check the enumeration using an independent F2 span oracle, all-vector stability/intersection tests, Frobenius pairing identities, the same-trace counterexample, the published torsion data and CLI reproducibility. The combined suite has 34 tests.

## Directions left open by this experiment

The most concrete next step is target reconstruction for additional verified kernels, followed by rationality, geometry and point-transfer checks. One should then test curve families not already covered by known attack analyses, rather than treating every positive kernel test as new. Higher-dimensional restrictions of genus-2, Picard or superelliptic Jacobians would require larger symplectic-module searches and an appropriate target reconstruction method; the current code does not perform those searches.

Prym representations and alternative polarizations are plausible search targets because they retain group structure while changing its geometric presentation. They remain hypotheses for a new speedup here. Arbitrary surfaces, moduli spaces or characteristic-zero lifts need an explicit compatible action or homomorphism before they become DLP reductions. No new general finite-field/number-field correspondence, practical attack, or NIST break has been established.

Follow-up computations are in [explicit quotient reconstruction](QUOTIENT_RECONSTRUCTION.md): 24 rational quartic constructions, 22 matching Frobenius polynomials, two opposite-sign subgroup obstructions, and actual toy divisor transfers. A source-equivalence audit identifies one underlying Weil-restriction case; the counts are constructions rather than independent families. Efficient divisor inversion remains open in this implementation.
