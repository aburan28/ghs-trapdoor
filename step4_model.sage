load("ghs_core.sage")
v, beta, bs, NE = load("weak_params.sobj")
set_random_seed(4)

K2 = GF(2^(2*N), 'c2')
fK = K.modulus()
roots = fK.change_ring(K2).roots()
img_a = roots[0][0]
emb = K.hom([img_a], K2)
def toK2(x): return emb(x)
def f2vec2(x):
    vv = list(x.polynomial()); return vector(F2, vv + [0]*(2*N-len(vv)))
_R2 = PolynomialRing(GF(2), 'X')
def f2from2(vv):
    return K2(_R2(list(vv)))
LM = Matrix(F2, [f2vec2(K2.gen()^(2*i) + K2.gen()^i) for i in range(2*N)]).transpose()
assert LM.rank() == 2*N - 1
def as_solve(c):
    return f2from2(LM.solve_right(f2vec2(c)))

# orbit data in K
vo = [v]
for i in range(3): vo.append(sigma(vo[-1]))
bo = [w^2 for w in vo]
sqrts = [sqrtK(bb) for bb in bo]        # v_i
# mu with Tr_{K/k}(mu) = 1
mu = None
for _ in range(200):
    cand = K.random_element()
    s = cand
    t = cand
    for i in range(3):
        t = sigma(t); s = s + t
    if s == K(1): mu = cand; break
assert mu is not None
mus = [mu]
for i in range(3): mus.append(sigma(mus[-1]))
# a = sqrt(gamma2/gamma1), gamma1 = v+v1, gamma2 = v+v2
g1 = vo[0]+vo[1]; g2 = vo[0]+vo[2]
a_cc = sqrtK(g2/g1)
acc = [a_cc]
for i in range(3): acc.append(sigma(acc[-1]))

def fiber(x0):
    """8 tuples (s0,s1,s2) with s_j^2+s_j = x0+v_j/x0, over K2; x0 given as K2 element"""
    sols = []
    for j in range(3):
        c = toK2(sqrts[j])/x0 + x0
        r = as_solve(c)
        sols.append([r, r+K2(1)])
    tuples = []
    for e0 in sols[0]:
        for e1 in sols[1]:
            for e2 in sols[2]:
                tuples.append((e0,e1,e2))
    return tuples

def invariants(tup):
    """given (s0,s1,s2) in K2, compute s3, ctilde, shat, u (all in K2)"""
    s0,s1,s2 = tup
    s3 = s0+s1+s2
    ss = [s0,s1,s2,s3]
    # c_i = (s_i + s_{i+2}) + a_i (s_i + s_{i+1})
    cs = []
    for i in range(4):
        ai = toK2(acc[i])
        cs.append(ss[i]+ss[(i+2)%4] + ai*(ss[i]+ss[(i+1)%4]))
    ctil = K2(0)
    shat = K2(0)
    for i in range(4):
        mi = toK2(mus[i])
        ctil += mi*cs[i]
        shat += mi*ss[i]
    u = shat^2 + shat
    return ctil, shat, u

# sanity: E-identity at a random tuple
x0 = toK2(k.random_element())
while x0 == 0: x0 = toK2(k.random_element())
tups = fiber(x0)
okY = True
for tup in tups[:4]:
    s0,s1,s2 = tup; s3 = s0+s1+s2
    w = [x0*s0+toK2(sqrts[0]), x0*s1+toK2(sqrts[1]), x0*s2+toK2(sqrts[2]), x0*s3+toK2(sqrts[3])]
    # y = w * T^{-1}: T[i][j] = th[(i+j)%4];  y_i = sum_j Tinv[j][i] w_j
    Tm = Matrix(K2, 4, 4, lambda i,j: toK2(th[(i+j)%4]))
    Tinv = Tm.inverse()
    ys = [sum(Tinv[j][i]*w[j] for j in range(4)) for i in range(4)]
    Y = sum(ys[i]*toK2(th[i]) for i in range(4))
    if Y^2 + x0*Y + x0^3 + toK2(beta) != 0: okY = False
print("E-identity on fiber points:", okY, " #tuples:", len(tups))

# collect (ctilde, u) pairs
pairs = []
xs_used = []
for t in range(120):
    xk = k.random_element()
    if xk == 0: continue
    x0 = toK2(xk)
    for tup in fiber(x0):
        ct, sh, uu = invariants(tup)
        pairs.append((ct, uu))
    xs_used.append(xk)
print("collected pairs:", len(pairs))
uniq = set()
for ct, uu in pairs: uniq.add((ct,uu))
print("unique (ctilde,u) pairs:", len(uniq))
save((pairs[:700], pairs[700:]), "fit_pairs.sobj")
save((v, beta, bs, NE, mu, acc), "params.sobj")
print("saved mu/acc into params.sobj")
