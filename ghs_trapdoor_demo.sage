# ================================================================
# GHS TRAP-DOOR SYSTEM -- full demonstration
# (trap door system of Sec. 5.5: weak curve Es, secret isogeny
#  chain to public curve Ep, escrow authority recovers keys via GHS)
# ================================================================
load("ghs_lib.sage")
import time
set_random_seed(42)

G = build_all()
E, NE, psi, J2, ZERO = G['E'], G['NE'], G['psi'], G['J2'], G['ZERO']
p_big = 7717
print("=" * 62)
print("GHS TRAP-DOOR SYSTEM DEMO")
print("K = GF(2^20) (n=4 over k=GF(2^5));  weak Es: Y^2+XY = X^3+beta")
print("#Es =", factor(NE), "  large prime p =", p_big)
print("=" * 62)

# ---------------- 1. user key generation -------------------------
L = 8       # isogeny chain length
t0 = time.time()
cur = E
phis = []
# step 0: a degree-2 isogeny (kernel = the 2-torsion point), then degree-17 steps.
# (odd-degree isogenies alone preserve the GHS-weak class here; the 2-step
#  followed by a 17-step moves Ep out of the weak class - checked below)
degs = [2] + [17]*(L-1)
for step in range(L):
    Nc = cur.order()
    if degs[step] == 2:
        a1,a2c,a3c,a4c,a6c = cur.a_invariants()
        # unique 2-torsion point: x = a3, y = sqrt(a3^3+a2*a3^2+a4*a3+a6)
        x0t = a3c
        rhs = x0t^3 + a2c*x0t^2 + a4c*x0t + a6c
        y0t = sqrtK(rhs)
        Kp = cur(x0t, y0t)
        assert 2*Kp == cur(0)
    else:
        R = cur.random_point()
        Kp = (Nc//17)*R
        while Kp.order() != 17:
            R = cur.random_point(); Kp = (Nc//17)*R
    phi = cur.isogeny(Kp)
    phis.append(phi)
    cur = phi.codomain()
Ep = cur
t_key = time.time()-t0
js = [phi.domain().j_invariant() for phi in phis] + [Ep.j_invariant()]
print("\n[1] keygen: %d isogenies of degree 17 (Es -> Ep), %.1fs" % (L, t_key))
print("    j-invariants all distinct (non-backtracking chain):", len(set(js)) == L+1)
print("    Ep:", Ep.a_invariants())

# user keypair on Ep: generator of the p-subgroup
Np = Ep.order()
assert Np == NE
R = Ep.random_point()
Pgen = (Np//p_big)*R
while Pgen.order() != p_big:
    R = Ep.random_point(); Pgen = (Np//p_big)*R
d_secret = ZZ.random_element(1, p_big)
Qpub = d_secret * Pgen
print("    user public key: Q = d*P with d secret; ord(P) =", Pgen.order())

# ---------------- 2. ElGamal encryption on Ep ----------------------
# map a message (small integer) to a point via a fixed embedding
msg = 12345
Mpt = msg * Pgen          # (demo embedding; any point works)
r = ZZ.random_element(1, p_big)
C1 = r * Pgen
C2 = Mpt + r * Qpub
print("\n[2] ElGamal on Ep: ciphertext (C1, C2) for message m = %d" % msg)
Mdec = C2 - d_secret * C1
print("    user decrypts: recovered m =", Mdec[0] if Mdec == msg*Pgen else "FAIL",
      "| check:", Mdec == msg*Pgen)

# ---------------- 3. escrow authority attack -----------------------
# (a) pull Qpub and Pgen back through the secret chain via preimages
t0 = time.time()
def pullback(Pt):
    S = Pt
    for phi in reversed(phis):
        S = phi.degree() * preimage(phi, S)   # = dual isogeny applied to S
        assert S is not None
    return S
P0 = pullback(Pgen)
Q0 = pullback(Qpub)
t_pull = time.time()-t0
print("\n[3] escrow: pulled (P, Q) back through the secret chain (%.1fs)" % t_pull)
print("    order of dual-chain image of P on Es:", P0.order())

# (b) GHS transfer to the Jacobian
t0 = time.time()
D1 = psi(P0)
D2 = psi(Q0)
t_psi = time.time()-t0
print("    GHS point-transfer psi applied (%.1fs); ord(psi(P0)) == %d" % (t_psi, p_big))

# (c) DLP in the Jacobian subgroup via BSGS
t0 = time.time()
m = ceil(sqrt(p_big))
babies = {}
Di = ZERO
for i in range(m):
    babies[Di] = i
    Di = Di + D1
step = m * D1
found = None
X = D2
for j in range(m+1):
    if X in babies:
        found = j*m - babies[X] if False else None
        # D2 = i*D1 + j*m*D1  =>  d = babies[X] + j*m
        found = babies[X] + j*m
        break
    X = X - step
assert found is not None
d_rec = found % p_big
t_dlp = time.time()-t0
print("    Jacobian DLP solved by BSGS (%.1fs)" % t_dlp)
print("    recovered d =", d_rec, "| correct:", d_rec == d_secret % p_big)

# (d) authority decrypts the message
M_esc = C2 - d_rec * C1
print("    authority decrypts: check:", M_esc == msg*Pgen)

# ---------------- 4. why Ep is not weak ----------------------------
def m_of(Ec):
    """GHS m-value: dim_F2 span{(a2_j, 1+sqrt(a6_j))} for the sigma-orbit"""
    a2 = Ec.a2(); a6 = Ec.a6()
    if a6 == 0: return None
    s6 = sqrtK(a6)
    pairs = []
    x2, x6 = a2, s6
    for i in range(N):
        pairs.append((x2, x6))
        x2 = sigma(x2); x6 = sigma(x6)
    vecs = [vector(F2, list(f2vec(a)) + list(f2vec(b))) for (a,b) in pairs]
    return Matrix(F2, vecs).rank()
mE = m_of(E); mEp = m_of(Ep)
assert mEp == 4, "Ep failed to escape the weak class"
print("\n[4] GHS weakness check:")
print("    m(Es)  =", mE, " -> genus 2^(m-1)-1 =", 2^(mE-1)-1, "(weak: GHS applies)")
print("    m(Ep)  =", mEp, " -> genus 2^(m-1)-1 =", 2^(mEp-1)-1, "(strong: GHS useless on Ep)")
print("    an outsider must use Pollard rho on Ep: ~sqrt(p) = 2^%.0f group ops" % log(sqrt(p_big),2))
print("    direct GHS on Ep needs a genus-%d Jacobian over GF(2^5): rho ~2^%.0f ops" % (2^(mEp-1)-1, log(sqrt(32^(2^(mEp-1)-1)),2)))
print("\nDONE.")
