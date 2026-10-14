# GHS Trap-Door System (key escrow via the GHS Weil-descent attack)

Implementation of the trap-door system described in the passage quoted in the task
("5.5. Trap Door Systems"): a user generates a secret elliptic curve **Es** that
is *susceptible to the GHS attack*, computes a public curve **Ep** by a secret,
long, random isogeny chain starting at Es, and submits Es together with the
chain to a trusted authority for key escrow. Ep is used as usual in elliptic
curve cryptosystems. Parameters are chosen so that Pollard-rho is the best
attack on Ep, while the DLP on Es is much easier via GHS Weil descent — so the
authority can decrypt with considerable but bounded computing effort.

Built on Gaudry–Hess–Smart, *Constructive and Destructive Facets of Weil
Descent on Elliptic Curves*, J. Cryptology 15(1), 2002 (paper included as
`refs_ghs_paper.pdf`), including the full conorm–norm point transfer
psi : E(K) -> Jac(H)(k) of its §3.6, implemented and verified from scratch in
SageMath 10.

## The construction (demo parameters)

* k = GF(2^5), K = GF(2^20) = k^4, sigma = Frobenius^5 (order 4), normal
  basis {theta, theta^2, theta^4, theta^3} with theta a primitive 5th root of
  unity (trace 1, squaring = cyclic shift — the "sparse" model of GHS §2).
* Weak curve **Es : Y^2 + XY = X^3 + beta**, beta = v^2, where v is sampled so
  that the GHS invariant m = dim_F2 Span{(1, sigma^i(v))} equals **3**.
  #Es = 2^3 * 17 * **7717** (7717 prime = the protected subgroup order p).
* GHS descent of Es: hyperelliptic curve **H / k : W^2 + Q(V) W = Pr(V) Q(V)**,
  genus **4 = 2^(m-1)** (verified by rigorous function-field genus computation).
* The point transfer psi (conorm from E to the AS-tower F, norm back to the
  fixed field F') is evaluated numerically as a sum of 16 conjugate points on
  H and reduced in the Jacobian. Verified: **homomorphism** (5+ random trials,
  psi(P1+P2) = psi(P1)+psi(P2)), and **order-preserving on the p-subgroup**
  (ord psi(P) = 7717), so the DLP transfers faithfully.
* Trap door: a secret chain Es -> Ep of one degree-2 and seven degree-17
  isogenies (all degrees coprime to p). Non-backtracking (distinct
  j-invariants). The chain is the escrowed secret.
* Escrow attack: pull (P, Q=dP) back through the chain using per-step
  preimages of the isogenies' x-rational maps and the identity
  phid(P) = deg(phi) * (any preimage) (dual isogenies are not implemented in
  char 2 in Sage — this avoids them entirely). Then transfer with psi and run
  baby-step/giant-step in the order-7717 subgroup of Jac(H)(GF(2^40)).

## Why the trap door closes

* m(Es) = 3 -> the GHS curve of Es has genus 4 over GF(2^5); the p-subgroup
  lands in Jac(H)(GF(2^5)) and the DLP there costs ~sqrt(p) — the authority's route.
* m(Ep) = 4 -> a direct GHS attack on Ep needs a genus-7 Jacobian over
  GF(2^5): rho there costs ~2^18 vs ~sqrt(p) = 2^6.5 with Pollard on Ep.
  An outsider uses Pollard on Ep; the authority uses GHS on Es.

An interesting empirical finding baked into the demo design: chains of
**odd**-degree (17-) isogenies alone *preserve* the weak class (m stayed 3 in
10/10 random chains — the trace conditions defining m=3 track the 2-Sylow
structure of E(K), which odd-degree isogenies preserve). A degree-2 step
followed by a degree-17 step escapes the class (m(Ep) = 4), which is why the
chain starts with a 2-isogeny. Note the related classical caveat
(Galbraith–Hess–Smart "Extending the GHS Weil descent attack"): an attacker who
can *find* an isogeny from Ep to a weak curve can extend the attack; at real
parameter scale (e.g. GF(2^155), as in the original trap-door proposal) the
weak class is sparse in a huge isogeny graph and the chain is long, which is
the intended countermeasure.

## Files

Generation pipeline (deterministic, run in order):
1. `step3a.sage` — search the weak curve (m = 3, Tr conditions, good group
   order) -> `weak_params.sobj`
2. `step4_model.sage` — field layer, Artin–Schreier solver, invariant
   evaluation on fibers; collects (ctilde, u) pairs -> `fit_pairs.sobj`,
   `params.sobj` (mu, a-coefficients, etc.)
3. `step5_fit.sage` — linear-algebra fit u = P(ctilde)/Q(ctilde) over GF(2^5)
   (extrapolation-validated) -> `H_model.sobj`
4. `step6_smooth.sage` — desingularize at infinity (Ep-reduction
   s -> s + alpha V^2 + beta V) -> smooth genus-4 model `H_smooth.sobj`

Library:
* `ghs_core.sage` — field setup, normal-basis coordinates, sigma, sqrt
* `ghs_lib.sage` — `build_all()` (curve, tower, psi), `preimage()` (isogeny
  preimage via x-rational map), AS solver

Demo and tests:
* `ghs_trapdoor_demo.sage` — the full trap-door system: keygen (isogeny
  chain), ElGamal on Ep, escrow (pullback + GHS + Jacobian DLP), weakness
  analysis. Output in `demo.txt`.
* `verify.sage` — verification battery (Weil-restriction identities,
  homomorphism property, order preservation, model genus). Output in
  `verify_out.txt`.

Run with SageMath (tested: Sage 10.10):

    sage ghs_trapdoor_demo.sage
    sage verify.sage

## Verification summary (all pass)

* sparse Weil-restriction equations y_{j-1}^2 + x y_j = x^3 + b_j (30 random)
* E-identity on all 16 fiber points of the AS-tower
* psi(P1+P2) = psi(P1) + psi(P2) in Jac(H) (random trials)
* ord psi(P) = 7717 = ord P on the protected subgroup
* genus(H) = 4 = 2^(m-1), matching the plane-component function field
* end-to-end: authority recovers d and decrypts the ElGamal ciphertext

## Scale and honest caveats

This is a *mechanism demo*: at GF(2^20) everything is tiny and Pollard on Ep
(~2^6.5 ops) and the authority's BSGS (~2^7 ops) are both trivial. The GHS
advantage is asymptotic: for real parameters (n = 155-type fields, as in the
trap-door literature) the authority's route is subexponential index calculus
on a moderate-genus Jacobian over GF(2^31) (Enge–Gaudry / Gaudry's
algorithm, O(q^2) for fixed genus) versus O(q^{n/2}) Pollard for outsiders.
The point-transfer machinery here (conorm–norm, order preservation, isogeny
pullback) is exactly the machinery that would be used at scale; only the DLP
solver and field sizes would change.

## Lean formalization (`lean/`)

A mathlib-free Lean 4 formalization of the core algebra is in `lean/`
(`GhsTrapdoor/Basic.lean`, builds with `cd lean && lake build`):

* `asp_add`, `asp_eq_zero`, `asp_fiber_shift` — the Artin–Schreier operator
  ℘(t) = t² + t over a characteristic-2 ring is additive, its kernel on a
  domain is {0,1}, and its fibres are the pairs {t, t+1} — the source of the
  two-component twist ε and the hyperelliptic involution ŝ ↦ ŝ+1 in the GHS
  construction.
* `trace_asp_zero` — Σ_{i<n} ℘(x)^{2^i} = 0 over any ring on which Frobenius
  x ↦ x^{2^n} is the identity: the trace-zero solvability condition behind
  the AS-solver.
* `tau_on`, `tau_tau` — the involution (x,y) ↦ (x, y+x) of
  E : y² + xy = x³ + β preserves the curve and is an involution: it
  generates the two components of the Weil restriction and the hyperelliptic
  structure of the descent curve.
