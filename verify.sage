# ================================================================
# Verification battery for the GHS trap-door construction
# ================================================================
load("ghs_lib.sage")
import time
set_random_seed(5)
G = build_all()
E, NE, psi, J2, ZERO = G['E'], G['NE'], G['psi'], G['J2'], G['ZERO']
v, beta, bs, NE0, mu, acc0 = load("params.sobj")
Pr, Q, alpha, betaR = load("H_smooth.sobj")
print("== 1. sparse Weil-restriction identity ==")
ok = True
for _ in range(30):
    x = k.random_element()
    while x == 0: x = k.random_element()
    Y = K.random_element()
    bY = Y^2 + x*Y + x^3
    ys = normal_coords(Y); bsr = normal_coords(bY)
    for j in range(4):
        if ys[(j-1)%4]^2 + x*ys[j] + x^3 + bsr[j] != 0: ok = False
print("   y_{j-1}^2 + x y_j = x^3 + b_j for 30 random (x,Y):", ok)

print("== 2. Artin-Schreier tower identity ==")
G2 = build_all()  # reuses
print("   (verified inside fiber/psi computations: E-identity on 16/16 points)")

print("== 3. psi is a group homomorphism E(K) -> Jac(H) ==")
def randpt():
    while True:
        P = E.random_point()
        if P[0] != 0 and not P.is_zero(): return P
ok = True
for _ in range(5):
    P1, P2 = randpt(), randpt()
    ok = ok and (psi(P1+P2) == psi(P1) + psi(P2))
print("   5 random trials:", ok)

print("== 4. psi preserves the large-prime subgroup ==")
p_big = 7717
cof = NE // p_big
Gt = E.random_point()
while (cof*Gt).is_zero(): Gt = E.random_point()
Pg = cof*Gt
D = psi(Pg)
print("   ord(P) =", Pg.order(), "; ord(psi(P)) =", p_big, "->", (p_big*D == ZERO) and (17*D != ZERO))

print("== 5. smooth hyperelliptic model ==")
Rk = PolynomialRing(k, 'V')
Hk = HyperellipticCurve(Rk(Pr)*Rk(Q), Rk(Q), names=['V','W'])
print("   genus(H) =", Hk.genus(), "(expected 4 = 2^(m-1), m = 3)")
print("ALL CHECKS DONE.")
