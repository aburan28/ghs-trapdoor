# Generation step 4/4: desingularize at infinity -> smooth genus-4 model
load("ghs_core.sage")
Pk, Qk = load("H_model.sobj")
Rk = PolynomialRing(k, 'V')
P, Q = Rk(Pk), Rk(Qk)
p = P.list(); q = Q.list()
def coeff(pl, i): return pl[i] if i < len(pl) else k(0)
sq = lambda c: c^(2^(l-1))
p8, p6, q4, q2 = coeff(p,8), coeff(p,6), coeff(q,4), coeff(q,2)
alpha = sq(p8/q4)
c6 = p6 + alpha^2*q2 + alpha*q4
betaR = sq(c6/q4)
Rr = alpha*Rk.gen()^2 + betaR*Rk.gen()
Pr = P + (Rr^2 + Rr)*Q
assert Pr.degree() == 5
Hk = HyperellipticCurve(Pr*Q, Q, names=['V','W'])
assert Hk.genus() == 4
print("smooth model: W^2 + Q(V)W = Pr(V)*Q(V), genus 4")
v, beta, bs, NE, mu, acc0 = load("params.sobj")
save((Pr, Q, alpha, betaR), "H_smooth.sobj")
print("saved H_smooth.sobj")
