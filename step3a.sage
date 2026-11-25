load("ghs_core.sage")
set_random_seed(3)
def orbit4(v):
    o = [v]
    for _ in range(3): o.append(sigma(o[-1]))
    return o
def m_dim_v(v):
    o = orbit4(v)
    vecs = [vector(F2, list(f2vec(K(1))) + list(f2vec(oo))) for oo in o]
    return Matrix(F2, vecs).rank()
def trace_Kk(v):
    return sum(orbit4(v))
import time
best = None
for _ in range(400):
    cand = K.random_element()
    if cand == 0 or trace_Kk(cand) != 0 or m_dim_v(cand) != 3: continue
    beta = cand^2
    E = EllipticCurve(K, [1,0,0,0,beta])
    NE = E.order()
    f = factor(NE)
    big = max([p for p,e in f])
    if best is None or big > best[0]:
        best = (big, cand, beta, NE, f)
    if big > 2^14: break
big, v, beta, NE, f = best
print("largest prime:", big, "=", factor(big))
print("#E =", f)
bs = normal_coords(beta)
print("bs:", bs)
save((v, beta, bs, NE), "weak_params.sobj")
