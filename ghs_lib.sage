# ============================================================
# GHS trap-door system library
# (Lenstra-style key-escrow trapdoor via the GHS Weil-descent
#  attack, Gaudry-Hess-Smart 2002; demo parameters n=4, l=5)
# ============================================================
load("ghs_core.sage")

def build_all():
    """returns dict with all GHS objects for the weak curve Es"""
    v, beta, bs, NE, mu, acc0 = load("params.sobj")
    Pr, Q, alpha, betaR = load("H_smooth.sobj")
    K2 = GF(2^(2*N), 'c2')
    fK = K.modulus()
    img_a = fK.change_ring(K2).roots()[0][0]
    emb = K.hom([img_a], K2)
    def E2(x): return emb(x)
    def k2(x): return emb(K(x))
    def f2vec2(x):
        vv = list(x.polynomial()); return vector(F2, vv + [0]*(2*N-len(vv)))
    _R2 = PolynomialRing(GF(2),'X')
    def f2from2(vv): return K2(_R2(list(vv)))
    LM = Matrix(F2, [f2vec2(K2.gen()^(2*i) + K2.gen()^i) for i in range(2*N)]).transpose()
    def as_solve(c): return f2from2(LM.solve_right(f2vec2(c)))
    vo = [v]
    for i in range(3): vo.append(sigma(vo[-1]))
    sqrts = vo
    mus = [mu]
    for i in range(3): mus.append(sigma(mus[-1]))
    acc = acc0
    Tm = Matrix(K2, 4, 4, lambda i,j: E2(th[(i+j)%4]))
    Tinv = Tm.inverse()
    K2V = PolynomialRing(K2, 'V')
    Q2 = K2V([k2(c) for c in Q.list()])
    Pr2 = K2V([k2(c) for c in Pr.list()])
    R2 = k2(alpha)*K2V.gen()^2 + k2(betaR)*K2V.gen()
    H2 = HyperellipticCurve(Pr2*Q2, Q2, names=['V','W'])
    J2 = H2.jacobian()
    ZERO = J2(0)
    E = EllipticCurve(K, [1,0,0,0,beta])

    def fiber_tuples(x0):
        sols = []
        for j in range(3):
            c = E2(sqrts[j])/x0 + x0
            r = as_solve(c)
            sols.append([r, r+K2(1)])
        return [(e0,e1,e2) for e0 in sols[0] for e1 in sols[1] for e2 in sols[2]]

    def Yval(x0, tup):
        s0,s1,s2 = tup; s3 = s0+s1+s2
        ss = [s0,s1,s2,s3]
        w = [x0*ss[j] + E2(sqrts[j]) for j in range(4)]
        ys = [sum(Tinv[j][i]*w[j] for j in range(4)) for i in range(4)]
        return sum(ys[i]*E2(th[i]) for i in range(4))

    def invariants(tup):
        s0,s1,s2 = tup; s3 = s0+s1+s2
        ss = [s0,s1,s2,s3]
        cs = []
        for i in range(4):
            ai = E2(acc[i])
            cs.append(ss[i]+ss[(i+2)%4] + ai*(ss[i]+ss[(i+1)%4]))
        ctil = K2(0); shat = K2(0)
        for i in range(4):
            mi = E2(mus[i])
            ctil += mi*cs[i]; shat += mi*ss[i]
        return ctil, shat

    def psi(P):
        """GHS point-transfer psi: E(K) -> Jac(H)(K2); verified homomorphism."""
        xP = E2(P[0]); yP = E2(P[1])
        tups = fiber_tuples(xP)
        sel = [t for t in tups if Yval(xP, t) == yP]
        assert len(sel) == 4, len(sel)
        D = ZERO
        for t in sel:
            ct, sh = invariants(t)
            Wr = Q2(ct)*(sh + R2(ct))
            ci, wi = ct, Wr
            for i in range(4):
                D = D + J2(H2([ci, wi]))
                ci = ci^(2^l); wi = wi^(2^l)
        return D

    return dict(E=E, NE=NE, psi=psi, J2=J2, ZERO=ZERO, H2=H2, K2=K2, emb=emb, beta=beta)

def preimage(phi, P):
    """find R in domain(phi)(K) with phi(R) = P (P must lie in the image)."""
    Edom = phi.domain()
    Kk = Edom.base_field()
    xP = P[0]
    xr = phi.x_rational_map()          # rational function in x
    Rr = xr.parent()
    X = Rr.gen()
    num = xr.numerator() - xP * xr.denominator()
    num = num/num.leading_coefficient()
    roots = num.roots(ring=Kk)
    for r, mlt in roots:
        # y from y^2 + x*y = x^3+a2x^2+a4x+a6
        a1,a2,a3,a4,a6 = Edom.a_invariants()
        c = r^3 + a2*r^2 + a4*r + a6
        # y^2 + r*y = c ; y = r*z ; r^2(z^2+z) = c
        if r == 0: continue
        cc = c / (r^2)
        tr = sum(cc^(2^i) for i in range(Kk.degree()))
        if tr != 0: continue
        z = as_solveK(Kk, cc)
        for zz in (z, z + Kk(1)):
            R = Edom(r, r*zz)
            if phi(R) == P:
                return R
    return None

_as_cache = {}
def as_solveK(Kk, c):
    """solve z^2+z=c in Kk (trace must be 0) via cached F2-linear algebra"""
    key = Kk
    if key not in _as_cache:
        F2 = GF(2)
        d = Kk.degree()
        def fv(x):
            vv = list(x.polynomial()); return vector(F2, vv + [0]*(d-len(vv)))
        RP = PolynomialRing(F2, 'X')
        M = Matrix(F2, [fv(Kk.gen()^(2*i) + Kk.gen()^i) for i in range(d)]).transpose()
        _as_cache[key] = (M, fv, d, RP)
    M, fv, d, RP = _as_cache[key]
    z = M.solve_right(fv(c))
    return Kk(RP(list(z)))
