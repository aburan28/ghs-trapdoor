"""Explicit theta quotients for the two verified F343 toy curves.

Requires python-flint==0.9.0 (research/requirements-reconstruction.txt).
All finite fields are exact. This module reconstructs plane quartics and
explicit quotient-to-Kummer maps. Jacobian divisor inversion is a separate
Sage experiment; no cryptographic speedup is asserted.
"""

import argparse
import itertools
import json
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx

from research.small_fields import frobenius_polynomial
from research.weil_quotient_kernels import (
    coordinates,
    frobenius_matrix,
    lagrangians,
    row_times,
)

if not __debug__:
    raise RuntimeError("The research certificates require Python assertions; omit -O")

# Explicit irreducible polynomial, checked by FLINT. F7^18 contains the
# F7^9 two-division splitting field and F7^6 three-torsion points.
MODULUS = [3, 2, 6, 0, 0, 3, 1, 5, 6, 1, 6, 2, 1, 0, 0, 0, 0, 0, 1]
K = fq_default_ctx(modulus=fmpz_mod_poly_ctx(7)(MODULUS), var="s")
R = fq_default_poly_ctx(K)
z = min(
    (a for a, m in R([1, 0, 1, 1]).roots()), key=lambda a: tuple(map(int, a.to_list()))
)


def encode(n):
    """Embed the census polynomial basis F7[z]/(z^3+z^2+1)."""
    if type(n) is not int or not 0 <= n < 343:
        raise ValueError("requires an integer field encoding in [0,343)")
    return K(n % 7) + K(n // 7 % 7) * z + K(n // 49) * z * z


def coefficients(value):
    return list(map(int, value.to_list()))


def canonical_sqrt(value):
    root = value.sqrt()
    return min((root, -root), key=lambda x: tuple(coefficients(x)))


def negative(P):
    return None if P is None else (P[0], -P[1])


def add(A, P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u:
        if y + v == 0:
            return None
        m = (3 * x * x + A) / (2 * y)
    else:
        m = (v - y) / (u - x)
    a = m * m - x - u
    return a, m * (x - a) - y


class EllipticWeilSet:
    """Exact normal 2/3-Weil functions; origin is shifted to a 2-torsion point.

    Raw 3-Weil functions are tangent lines at -P, evaluated at R-Q3.
    Normalization uses alpha_P=-1/c(P,P), verified against every pair of
    E[3] points. No Riemann-Roch black box is used in this stage.
    """

    def __init__(self, A, B, x1, x2):
        self.A, self.B = A, B
        roots = sorted(
            (x for x, m in R([B, A, 0, 1]).roots()),
            key=lambda a: tuple(map(int, a.to_list())),
        )
        assert len(roots) == 3
        self.two = [None] + [(x, K(0)) for x in roots]
        U, V = (
            (x1, canonical_sqrt(x1**3 + A * x1 + B)),
            (x2, canonical_sqrt(x2**3 + A * x2 + B)),
        )
        assert add(A, add(A, U, U), U) is None
        assert add(A, add(A, V, V), V) is None
        self.three = {}
        for a, b in itertools.product(range(3), repeat=2):
            aa = None if a == 0 else U if a == 1 else negative(U)
            bb = None if b == 0 else V if b == 1 else negative(V)
            self.three[a, b] = add(A, aa, bb)
        assert len(set(self.three.values())) == 9
        for n in range(50):
            x = K.gen() + n
            yy = x**3 + A * x + B
            if yy.is_square():
                self.probe = (x, canonical_sqrt(yy))
                break
        self.alpha = {(0, 0): K(1)}
        for ab, P in self.three.items():
            if P is None:
                continue
            pab = tuple(-a % 3 for a in ab)
            c = (
                self.raw(ab, self.probe)
                * self.raw(ab, add(A, self.probe, P))
                / self.raw(pab, self.probe)
            )
            self.alpha[ab] = -1 / c
        U, V = self.three[1, 0], self.three[0, 1]
        self.zeta = (
            self.w3((1, 0), self.probe)
            * self.w3((0, 1), add(A, self.probe, U))
            / (self.w3((0, 1), self.probe) * self.w3((1, 0), add(A, self.probe, V)))
        )
        assert self.zeta**3 == 1 and self.zeta != 1
        for ab, cd in itertools.product(self.three, repeat=2):
            pq = tuple((a + b) % 3 for a, b in zip(ab, cd))
            c = (
                self.w3(ab, self.probe)
                * self.w3(cd, add(A, self.probe, self.three[ab]))
                / self.w3(pq, self.probe)
            )
            assert c == self.zeta ** (2 * (ab[0] * cd[1] - ab[1] * cd[0]) % 3), (
                "3normal",
                ab,
                cd,
                c,
            )
        e1, e2, e3 = roots
        self.g = [
            K(1),
            canonical_sqrt((e1 - e3) / (e1 - e2)),
            canonical_sqrt((e2 - e3) / (e2 - e1)),
            canonical_sqrt((e3 - e1) * (e2 - e3)),
        ]
        ratio = (
            self.w2(1, self.probe)
            * self.w2(2, add(A, self.probe, self.two[1]))
            / self.w2(3, self.probe)
        )
        assert ratio in [K(1), K(-1)]
        if ratio == 1:
            self.g[3] = -self.g[3]
        for i, j in itertools.product(range(4), repeat=2):
            ratio = (
                self.w2(i, self.probe)
                * self.w2(j, add(A, self.probe, self.two[i]))
                / self.w2(i ^ j, self.probe)
            )
            expected = -1 if (i & 1) and (j & 2) else 1
            assert ratio == expected, ("2normal", i, j, ratio, expected)

    def raw(self, ab, point):
        P = self.three[ab]
        if P is None:
            return K(1)
        S = add(self.A, point, self.two[3])
        x, y = negative(P)
        assert S is not None
        m = (3 * x * x + self.A) / (2 * y)
        return S[1] - y - m * (S[0] - x)

    def w3(self, ab, point):
        return self.alpha[ab] * self.raw(ab, point)

    def w2(self, i, point):
        if i == 0:
            return K(1)
        if point is None:
            return self.g[i] if i != 3 else K(0)
        x = point[0]
        e1, e2, e3 = [P[0] for P in self.two[1:]]
        if i == 1:
            return self.g[i] * (x - e2) / (x - e3)
        if i == 2:
            return self.g[i] * (x - e1) / (x - e3)
        return self.g[3] / (x - e3)

    def conjugate(self, n):
        out = object.__new__(EllipticWeilSet)
        frob = lambda a: a.frobenius(n)
        point = lambda p: None if p is None else tuple(frob(x) for x in p)
        out.A, out.B = frob(self.A), frob(self.B)
        out.two = list(map(point, self.two))
        out.three = {ab: point(p) for ab, p in self.three.items()}
        out.alpha = {ab: frob(a) for ab, a in self.alpha.items()}
        out.g = list(map(frob, self.g))
        out.zeta = frob(self.zeta)
        return out


def span(rows):
    return [
        tuple(sum(c[i] * rows[i][j] for i in range(3)) % 3 for j in range(6))
        for c in itertools.product(range(3), repeat=3)
    ]


def theta_data(weils, rows, seed=(0, 0, 0), with_gradient=False):
    """Theta-null squares and odd-theta gradients in the source Lie coordinates.

    Project sum(f_P * translated(seed)) onto the kernel-invariant section.
    Seeds are tensor products of 1,x,y in L(3Q3). Constant-seed projection
    can vanish identically, which is not a geometric obstruction. At Q3,
    use the pole-cancelling first jet in the parameter t=-x/y at infinity.
    """
    elements = span(rows)
    values, derivatives = {}, {}
    for i, w in enumerate(weils):
        for ti in range(4):
            for ab, P in w.three.items():
                ss = seed[i]
                if ti == 3:
                    if P is None:
                        value, deriv = (
                            (K(-1), K(0))
                            if ss == 2
                            else (K(0), K(1) if ss == 1 else K(0))
                        )
                    else:
                        x, y = P
                        sp = K(1) if ss == 0 else x if ss == 1 else y
                        dp = K(0) if ss == 0 else 2 * y if ss == 1 else 3 * x * x + w.A
                        minus_slope = -(3 * x * x + w.A) / (2 * y)
                        value, deriv = (
                            -w.alpha[ab] * sp,
                            w.alpha[ab] * (-minus_slope * sp - dp),
                        )
                else:
                    S = add(w.A, w.two[ti], w.two[3])
                    U = add(w.A, S, P)
                    if P is None:
                        fp, dfp = K(1), K(0)
                    else:
                        fp, dfp = (
                            w.w3(ab, w.two[ti]),
                            w.alpha[ab] * (3 * S[0] * S[0] + w.A),
                        )
                    sp = K(1) if ss == 0 else U[0] if ss == 1 else U[1]
                    dp = (
                        K(0)
                        if ss == 0
                        else 2 * U[1]
                        if ss == 1
                        else 3 * U[0] * U[0] + w.A
                    )
                    value, deriv = fp * sp, dfp * sp + fp * dp
                values[i, ti, ab], derivatives[i, ti, ab] = value, deriv

    def lead(T):
        total = K(0)
        grad = [K(0)] * 3
        for P in elements:
            val = [values[i, T[i], P[2 * i : 2 * i + 2]] for i in range(3)]
            der = [derivatives[i, T[i], P[2 * i : 2 * i + 2]] for i in range(3)]
            total += val[0] * val[1] * val[2]
            for i in range(3):
                grad[i] += der[i] * val[(i + 1) % 3] * val[(i + 2) % 3]
        return total, grad

    out = []
    gradients = {}
    for n in range(64):
        T = [((n >> (3 + i)) & 1) + 2 * ((n >> i) & 1) for i in range(3)]
        g = K(1)
        for i, w in enumerate(weils):
            g *= w.g[T[i]]
        value, grad = lead(T)
        out.append(g**3 * value**2)
        if value == 0 and any(c != 0 for c in grad):
            gradients[n] = grad
    return (out, gradients) if with_gradient else out


def solve(a, b):
    a = [list(row) + [v] for row, v in zip(a, b)]
    n = len(a)
    for j in range(n):
        i = next(i for i in range(j, n) if a[i][j] != 0)
        a[j], a[i] = a[i], a[j]
        d = a[j][j]
        a[j] = [x / d for x in a[j]]
        for i in range(n):
            if i == j:
                continue
            c = a[i][j]
            a[i] = [x - c * y for x, y in zip(a[i], a[j])]
    return [row[-1] for row in a]


def poly_add(a, b, scale=1):
    out = dict(a)
    for mon, c in b.items():
        out[mon] = out.get(mon, K(0)) + scale * c
    return {mon: c for mon, c in out.items() if c != 0}


def poly_mul(a, b):
    out = {}
    for p, c in a.items():
        for q, d in b.items():
            r = tuple(x + y for x, y in zip(p, q))
            out[r] = out.get(r, K(0)) + c * d
    return {mon: c for mon, c in out.items() if c != 0}


def linear(row):
    return {
        tuple(int(i == j) for j in range(3)): c for i, c in enumerate(row) if c != 0
    }


def quartic(th):
    def t(*args):
        value = K(1)
        for i in args:
            value *= th[i]
        return value

    a = [
        K(1),
        (t(5, 12, 33, 40) + t(21, 28, 49, 56) - t(7, 14, 35, 42))
        / (2 * t(5, 12, 49, 56)),
        (t(5, 12, 33, 40) + t(7, 14, 35, 42) - t(21, 28, 49, 56))
        / (2 * t(5, 12, 35, 42)),
    ]
    b = [
        K(1),
        (t(2, 28, 47, 49) + t(5, 27, 40, 54) - t(14, 16, 35, 61))
        / (2 * t(5, 27, 47, 49)),
        (t(2, 28, 47, 49) - t(5, 27, 40, 54) - t(14, 16, 35, 61))
        / (2 * t(5, 27, 35, 61)),
    ]
    c = [
        K(1),
        (t(7, 16, 42, 61) - t(12, 27, 33, 54) - t(2, 21, 47, 56))
        / (-2 * t(12, 27, 47, 56)),
        (t(12, 27, 33, 54) + t(7, 16, 42, 61) - t(2, 21, 47, 56))
        / (-2 * t(12, 27, 42, 61)),
    ]
    al = [a, b, c]
    return riemann(al)


def riemann(al):
    """Riemann's quartic from a normalized Aronhold system of seven lines."""
    u = solve([[1 / al[j][i] for j in range(3)] for i in range(3)], [K(-1)] * 3)
    kk = solve([[u[j] * al[j][i] for j in range(3)] for i in range(3)], [K(-1)] * 3)
    lhs = [[K(1)] * 3] + [[1 / x for x in row] for row in al]
    rhs = [[K(-1)] * 3] + [[-kk[i] * x for x in row] for i, row in enumerate(al)]
    xi = [solve(lhs[:3], [rhs[j][i] for j in range(3)]) for i in range(3)]
    xi = list(map(list, zip(*xi)))
    assert all(
        sum(lhs[3][j] * xi[j][i] for j in range(3)) == rhs[3][i] for i in range(3)
    ), "Riemann system inconsistent"
    xx = [linear([K(int(i == j)) for j in range(3)]) for i in range(3)]
    zz = [linear(row) for row in xi]
    terms = [poly_mul(x, y) for x, y in zip(xx, zz)]
    s = poly_add(poly_add(terms[0], terms[1]), terms[2], -1)
    f = poly_add(poly_mul(s, s), poly_mul(terms[0], terms[1]), -4)
    return f, al, xi


def substitute(f, forms):
    out = {}
    for mon, c in f.items():
        term = {(0, 0, 0): c}
        for n, form in zip(mon, forms):
            for _ in range(n):
                term = poly_mul(term, linear(form))
        out = poly_add(out, term)
    return out


def gradient_quartic(gradients, indices=(63, 25, 44, 11, 22, 37, 50)):
    """Reconstruct directly in a rational basis of Lie(Res E).

    The fixed odd characteristics form an Aronhold system. This uses
    semilinear Galois descent of the Lie basis, not the differential of
    the inseparable Frobenius endomorphism.
    """
    lines = [gradients[i] for i in indices]
    basis = lines[:3]
    transpose = list(map(list, zip(*basis)))
    scale = solve(transpose, lines[3])
    basis = [[scale[i] * x for x in row] for i, row in enumerate(basis)]
    transpose = list(map(list, zip(*basis)))
    al = []
    for line in lines[4:]:
        c = solve(transpose, line)
        al.append([x / c[0] for x in c])
    f, al, _xi = riemann(al)
    # Source Weil-restriction Lie basis: (a,a^7,a^49), a=c0+c1*z+c2*z^2.
    base = [
        [sum(row[i] * z ** (j * 7**i) for i in range(3)) for j in range(3)]
        for row in basis
    ]
    out = substitute(f, base)
    lead = out[min(out)]
    out = {m: c / lead for m, c in out.items()}
    return out


def determinant(matrix):
    a = [list(row) for row in matrix]
    n = len(a)
    prev = R.one()
    sign = 1
    for j in range(n - 1):
        pivot = next((i for i in range(j, n) if a[i][j] != 0), None)
        if pivot is None:
            return R.zero()
        if pivot != j:
            a[pivot], a[j] = a[j], a[pivot]
            sign = -sign
        value = a[j][j]
        for i in range(j + 1, n):
            for k in range(j + 1, n):
                numerator = a[i][k] * value - a[i][j] * a[j][k]
                q, r = divmod(numerator, prev)
                assert r == 0
                a[i][k] = q
        for i in range(j + 1, n):
            a[i][j] = R.zero()
        prev = value
    return sign * a[-1][-1]


def resultant(a, b):
    while a and a[-1] == 0:
        a = a[:-1]
    while b and b[-1] == 0:
        b = b[:-1]
    m, n = len(a) - 1, len(b) - 1
    matrix = []
    for i in range(n):
        matrix.append([R.zero()] * i + list(reversed(a)) + [R.zero()] * (n - i - 1))
    for i in range(m):
        matrix.append([R.zero()] * i + list(reversed(b)) + [R.zero()] * (m - i - 1))
    return determinant(matrix)


def derivative(f, j):
    return {
        tuple(x - (int(i == j)) for i, x in enumerate(mon)): c * mon[j]
        for mon, c in f.items()
        if mon[j]
    }


def smooth(f):
    """Exact affine resultant and infinity certificates; False is inconclusive.

    Possible common resultant roots are checked by a third polynomial gcd,
    since coincidences between different fibers can be extraneous.
    """

    def as_y(poly):
        a = [R.zero() for _ in range(5)]
        for (i, j, k), c in poly.items():
            a[j] += R([K(0)] * i + [c])
        return a

    fx, fy, fz = [derivative(f, i) for i in range(3)]
    a, b = resultant(as_y(f), as_y(fy)), resultant(as_y(fx), as_y(fy))
    gg = a.gcd(b)
    if gg.degree() > 0:
        roots = gg.roots()
        if sum(n for x, n in roots) != gg.degree():
            return False
        for x, n in roots:
            common = R.zero()
            for poly in [f, fx, fy]:
                value = R.zero()
                for (i, j, k), c in poly.items():
                    value += R([K(0)] * j + [c * x**i])
                common = common.gcd(value)
            if common.degree() > 0:
                return False
    infinity = R.zero()
    for df in [fx, fy, fz]:
        p = R.zero()
        for (i, j, k), c in df.items():
            if k == 0:
                p += R([K(0)] * j + [c])
        infinity = infinity.gcd(p)
    origin = any(df.get((0, 3, 0), K(0)) != 0 for df in [fx, fy, fz])
    return infinity.degree() == 0 and origin


def restrict_line(f, line):
    pivot = next(i for i, x in enumerate(line) if x != 0)
    others = [i for i in range(3) if i != pivot]
    substitutions = [None] * 3
    substitutions[others[0]] = R([0, 1])
    substitutions[others[1]] = R.one()
    substitutions[pivot] = (
        -(line[others[0]] * substitutions[others[0]] + line[others[1]]) / line[pivot]
    )
    out = R.zero()
    for mon, c in f.items():
        term = R([c])
        for i, n in enumerate(mon):
            term *= substitutions[i] ** n
        out += term
    return out


@lru_cache(maxsize=3)
def counting_field(degree):
    from research.small_fields import SmallField

    k = SmallField(7, degree)
    k.tables()
    return k


def count_projective_prime(f, degree):
    k = counting_field(degree)
    add, mul = k.tables()
    size = k.size
    ff = {m: int(c.to_list()[0]) for m, c in f.items()}
    assert all(c.frobenius(1) == c for c in f.values())
    count = 0
    for x in range(size):
        powers = [1]
        for _ in range(4):
            powers.append(mul[powers[-1]][x])
        c = [0] * 5
        for (i, j, l), a in ff.items():
            c[j] = add[c[j]][mul[a][powers[i]]]
        for y in range(size):
            value = c[4]
            for j in [3, 2, 1, 0]:
                value = add[mul[value][y]][c[j]]
            count += value == 0
        value = 0
        for (i, j, l), a in ff.items():
            if l == 0:
                value = add[value][mul[a][powers[i]]]
        count += value == 0
    return count + int(ff.get((4, 0, 0), 0) == 0)


def kummer_sections(weils, rows, seed, point):
    """64 homogeneous quotient-to-Kummer sections g_T^3 * theta(R+T)^2.

    They avoid division by theta(R), including when it vanishes. Coordinates
    are over F7^18. Exceptional poles away from the evaluated subgroup and
    its 3-torsion translates are outside this bounded implementation.
    """
    elements = span(rows)
    values, gvalues = {}, {}
    for i, w in enumerate(weils):
        for ti in range(4):
            if point[i] is None:
                gvalues[i, ti] = w.g[ti]
            else:
                gvalues[i, ti] = w.w2(ti, point[i])
            for ab, P in w.three.items():
                ss = seed[i]
                if point[i] is None and ti == 3:
                    if P is None:
                        value = K(-1) if ss == 2 else K(0)
                    else:
                        value = -w.alpha[ab] * (
                            K(1) if ss == 0 else P[0] if ss == 1 else P[1]
                        )
                else:
                    shifted = add(w.A, point[i], w.two[ti])
                    S = add(w.A, shifted, w.two[3])
                    U = add(w.A, S, P)
                    sp = K(1) if ss == 0 else U[0] if ss == 1 else U[1]
                    value = w.w3(ab, shifted) * sp
                values[i, ti, ab] = value
    out = []
    for n in range(64):
        T = [((n >> (3 + i)) & 1) + 2 * ((n >> i) & 1) for i in range(3)]
        total = K(0)
        for P in elements:
            total += (
                values[0, T[0], P[:2]]
                * values[1, T[1], P[2:4]]
                * values[2, T[2], P[4:]]
            )
        out.append(
            (gvalues[0, T[0]] * gvalues[1, T[1]] * gvalues[2, T[2]]) ** 3
            * total
            * total
        )
    return out


def normalize(v):
    lead = next(x for x in v if x != 0)
    return tuple(x / lead for x in v)


def rational_bitangents(gradients):
    return {
        n: [sum(row[i] * z ** (j * 7**i) for i in range(3)) for j in range(3)]
        for n, row in gradients.items()
    }


def choose_seed(weils, rows):
    attempts = []
    for seed in itertools.product(range(3), repeat=3):
        if sum(x == 2 for x in seed) % 2:
            continue  # Only even sections can project to the required even theta.
        theta, gradients = theta_data(weils, rows, seed, True)
        nonzero = sum(x != 0 for x in theta)
        attempts.append({"seed": seed, "nonzero_theta_nulls": nonzero})
        if nonzero:
            assert nonzero == 36 and len(gradients) == 28
            return seed, theta, gradients, attempts
    raise ArithmeticError("No nonzero invariant section among the 14 even tensor seeds")


def subgroup_map_report(weils, rows, seed, generator_y):
    """Evaluate all 373 source points, then check all kernel translates at 3 points."""
    generator = (K(0), encode(generator_y))
    assert generator[1] ** 2 == weils[0].B
    point, images, digest = None, [], sha256()
    samples, selected = [], (0, 1, 2, 4, 8, 16, 32)
    for n in range(373):
        triple = tuple(
            None if point is None else tuple(x.frobenius(i) for x in point)
            for i in range(3)
        )
        image = normalize(kummer_sections(weils, rows, seed, triple))
        images.append(image)
        digest.update(bytes(c for x in image for c in coefficients(x)))
        if n in (0, 1, 2, 17, 37, 186, 372):
            small = normalize([image[i] for i in selected])
            samples.append(
                {
                    "source_scalar_for_validation": n,
                    "seven_coordinate_indices": selected,
                    "seven_projective_coordinates": [coefficients(x) for x in small],
                }
            )
        point = add(weils[0].A, point, generator)
    assert point is None and len(set(images)) == 187
    assert all(images[n] == images[-n % 373] for n in range(373))
    for n in (1, 17, 186):
        point = None
        for _ in range(n):
            point = add(weils[0].A, point, generator)
        triple = tuple(tuple(x.frobenius(i) for x in point) for i in range(3))
        for vector in span(rows):
            shifted = tuple(
                add(w.A, triple[i], w.three[vector[2 * i : 2 * i + 2]])
                for i, w in enumerate(weils)
            )
            assert normalize(kummer_sections(weils, rows, seed, shifted)) == images[n]
    return {
        "codomain": "Kummer quotient, coordinates over F7^18",
        "source_points_evaluated": 373,
        "distinct_projective_images": 187,
        "only_sign_pairs_collide_verified": True,
        "kernel_translation_checks": 81,
        "ordered_64_coordinate_images_sha256": digest.hexdigest(),
        "jacobian_divisor_lift_verified": False,
        "samples": samples,
    }


def run(include_subgroup_maps=True):
    census_path = Path(__file__).with_name("results") / "cubic-curve-census.json"
    census = json.loads(census_path.read_text())
    sources = next(f for f in census["fields"] if f["p"] == 7)[
        "prime_order_extension_j_candidates"
    ]
    frob = frobenius_matrix(7, 3, [[2, 0], [0, 2]])
    sources_out = []
    for source in sources:
        x1, x2 = source["three_division_roots"][:2]
        w = EllipticWeilSet(
            encode(source["a"]), encode(source["b"]), encode(x1), encode(x2)
        )
        weils = [w, w.conjugate(1), w.conjugate(2)]
        targets = []
        map_data = None
        for index, (rows, intersections) in enumerate(lagrangians(3, 1)):
            if any(intersections) or any(
                coordinates(row_times(row, frob, 3), rows, 3) is None for row in rows
            ):
                continue
            seed, _theta, gradients, attempts = choose_seed(weils, rows)
            curve = gradient_quartic(gradients)
            assert all(c.frobenius(1) == c for c in curve.values())
            assert smooth(curve), "Smoothness certificate inconclusive"
            bitangents = rational_bitangents(gradients)
            for line in bitangents.values():
                value = restrict_line(curve, line)
                scalar, factors = value.factor()
                assert all(n % 2 == 0 for factor, n in factors)
                assert (4 - value.degree()) % 2 == 0 and scalar.is_square()
            counts = [count_projective_prime(curve, n) for n in (1, 2, 3)]
            polynomial = frobenius_polynomial(7, counts)
            matches = polynomial == [343, 0, 0, 29, 0, 0, 1]
            assert matches or polynomial == [343, 0, 0, -29, 0, 0, 1]
            targets.append(
                {
                    "kernel_enumeration_index": index,
                    "kernel_basis_rows": rows,
                    "theta_seed": seed,
                    "seed_attempts": attempts,
                    "nonzero_theta_nulls": 36,
                    "bitangent_lines_checked": 28,
                    "quartic_coefficients_ascending_monomials": [
                        [m, int(c.to_list()[0])] for m, c in sorted(curve.items())
                    ],
                    "field": "F7",
                    "smoothness_certified": True,
                    "point_counts_over_F7_F49_F343": counts,
                    "frobenius_polynomial_ascending": polynomial,
                    "jacobian_order": sum(polynomial),
                    "matches_Weil_restriction_Frobenius": matches,
                    "status": "matching_isogeny_class_requires_divisor_map"
                    if matches
                    else "opposite_sign_Frobenius_rational_373_subgroup_obstructed",
                }
            )
            if map_data is None and matches:
                map_data = (rows, seed, source["subgroup_witness"]["point"][1], index)
        assert len(targets) == 12
        report = {
            "source_a": source["a"],
            "source_b": source["b"],
            "source_order": 373,
            "source_torsion_x_coordinates": [x1, x2],
            "targets": targets,
            "matching_targets": sum(
                t["matches_Weil_restriction_Frobenius"] for t in targets
            ),
            "opposite_sign_targets": sum(
                not t["matches_Weil_restriction_Frobenius"] for t in targets
            ),
        }
        if include_subgroup_maps:
            rows, seed, generator_y, index = map_data
            report["selected_kernel_subgroup_map"] = {
                "kernel_enumeration_index": index,
                **subgroup_map_report(weils, rows, seed, generator_y),
            }
        sources_out.append(report)
    return {
        "scope": {
            "toy_only": True,
            "known_theta_quotient_method": True,
            "explicit_curves_over_F7_constructed": True,
            "source_torsion_and_normal_Weil_sets_verified": True,
            "theta_derivative_reconstruction_implemented": True,
            "Jacobian_divisor_transfer_implemented_in_this_module": False,
            "practical_DLP_speedup_established": False,
            "novel_vulnerability_established": False,
        },
        "arithmetic": {
            "dependency": "python-flint==0.9.0",
            "characteristic": 7,
            "working_extension_degree": 18,
            "modulus_ascending": MODULUS,
            "census_z_embedding_coefficients": coefficients(z),
        },
        "sources": sources_out,
    }


def export_divisor_inputs():
    """Exact inputs for an independent Sage divisor lift on one target per source.

    Image vectors come from evaluating the quotient sections, without using
    any Jacobian DLP or the scalar labels used for subsequent validation.
    Hex strings encode 7*18 coefficients in [0,7), in the stated order.
    """
    sources = []
    reverse = {encode(n): n for n in range(343)}
    for a, b, x1, x2, generator_y in ((248, 42, 19, 65, 176), (3, 76, 36, 47, 119)):
        w = EllipticWeilSet(encode(a), encode(b), encode(x1), encode(x2))
        weils = [w, w.conjugate(1), w.conjugate(2)]
        rows = lagrangians(3, 1)[153][0]
        seed, theta, gradients, _ = choose_seed(weils, rows)
        curve = gradient_quartic(gradients)
        point = None
        generator = (K(0), encode(generator_y))
        vectors, source_points = [], []
        for _ in range(373):
            triple = tuple(
                None if point is None else tuple(x.frobenius(i) for x in point)
                for i in range(3)
            )
            values = kummer_sections(weils, rows, seed, triple)
            small = normalize([values[i] for i in (0, 1, 2, 4, 8, 16, 32)])
            vectors.append(bytes(c for x in small for c in coefficients(x)).hex())
            source_points.append(None if point is None else [reverse[c] for c in point])
            point = add(w.A, point, generator)
        assert point is None and len(set(vectors)) == 187
        sources.append(
            {
                "source_a": a,
                "source_b": b,
                "source_generator": [0, generator_y],
                "source_order": 373,
                "kernel_enumeration_index": 153,
                "kernel_basis_rows": rows,
                "theta_seed": seed,
                "quartic": [[m, int(c.to_list()[0])] for m, c in sorted(curve.items())],
                "bitangents": {
                    str(n): [coefficients(x) for x in line]
                    for n, line in rational_bitangents(gradients).items()
                },
                "theta_nulls": [coefficients(x / theta[0]) for x in theta],
                "source_points_for_validation": source_points,
                "source_seven_coordinate_images_hex": vectors,
            }
        )
    return {
        "working_field_modulus_ascending": MODULUS,
        "characteristic": 7,
        "working_extension_degree": 18,
        "coefficient_encoding": "7 elements, each 18 base-7 coefficients, hex bytes",
        "image_coordinate_indices": [0, 1, 2, 4, 8, 16, 32],
        "census_z_embedding_coefficients": coefficients(z),
        "sources": sources,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--export-divisor-inputs",
        action="store_true",
        help="Export quotient evaluations for the separate Sage divisor experiment",
    )
    parser.add_argument(
        "--skip-subgroup-maps",
        action="store_true",
        help="Reconstruct all targets without the two exhaustive subgroup evaluations",
    )
    args = parser.parse_args()
    data = (
        export_divisor_inputs()
        if args.export_divisor_inputs
        else run(not args.skip_subgroup_maps)
    )
    result = json.dumps(data, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
