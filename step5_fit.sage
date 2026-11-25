# Generation step 3/4: fit u = P(ctilde)/Q(ctilde) from the collected pairs
# (deterministic: only needs fit_pairs.sobj)
load("ghs_core.sage")
K2 = GF(2^(2*N), 'c2')
train, test = load("fit_pairs.sobj")
R.<V> = PolynomialRing(K2)
rows = [[ct^i for i in range(21)] + [uu*ct^j for j in range(13)] for (ct,uu) in train]
ns = Matrix(K2, rows).right_kernel().basis()
vec = ns[0]
P0 = R(list(vec[:21])); Q0 = R(list(vec[21:]))
g = P0.gcd(Q0); Pn = P0//g; Qn = Q0//g
Qn = Qn/Qn.leading_coefficient(); Pn = Pn/Qn.leading_coefficient()
assert all(Pn(ct) + uu*Qn(ct) == 0 for (ct,uu) in test)
def ink(x): return x^(2^l) == x
assert all(ink(c) for c in Pn.coefficients()) and all(ink(c) for c in Qn.coefficients())
print("fit OK: deg P =", Pn.degree(), ", deg Q =", Qn.degree())
# to k-polynomials
fK = K.modulus()
img_a = fK.change_ring(K2).roots()[0][0]
emb = K.hom([img_a], K2)
lut = {emb(K(e)): e for e in k}
_Rf2 = PolynomialRing(GF(2), 'Z')
def tok(c): return k(_Rf2(lut[c].polynomial().list()))
Rk = PolynomialRing(k, 'V')
Pk = Rk([tok(c) for c in Pn.list()]); Qk = Rk([tok(c) for c in Qn.list()])
save((Pk, Qk), "H_model.sobj")
print("saved H_model.sobj")
