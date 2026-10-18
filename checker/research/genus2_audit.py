"""Bounded branch-locus audit; no DLP or general Weil-descent implementation.

Run from the project root: python research/genus2_audit.py --output PATH
Only the standard library is required. Integers encode a+b*i as a+p*b;
the auxiliary splitting field is K[v]/(v^2-delta), encoded a+|K|*b.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path


def prime64(n):
    """Deterministic Miller-Rabin for 0 <= n < 2^64."""
    if not 0 <= n < 2**64:
        raise ValueError("prime must be smaller than 2^64")
    for d in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % d == 0:
            return n == d
    if n < 2:
        return False
    s, d = 0, n - 1
    while d % 2 == 0:
        s, d = s + 1, d // 2
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


class Fp2:
    def __init__(self, p):
        if not prime64(p) or p % 4 != 3:
            raise ValueError("requires prime p = 3 mod 4")
        self.p, self.size = p, p * p

    def elt(self, a, b=0):
        return a % self.p + self.p * (b % self.p)

    def add(self, a, b):
        p = self.p
        return (a % p + b % p) % p + p * ((a // p + b // p) % p)

    def neg(self, a):
        return self.elt(-a % self.p, -(a // self.p))

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        p = self.p
        ar, ai, br, bi = a % p, a // p, b % p, b // p
        return (ar * br - ai * bi) % p + p * ((ar * bi + ai * br) % p)

    def inv(self, a):
        if not a:
            raise ZeroDivisionError
        p = self.p
        ar, ai = a % p, a // p
        scale = pow((ar * ar + ai * ai) % p, -1, p)
        return self.elt(ar * scale, -ai * scale)

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def pow(self, a, n):
        if n < 0:
            return self.pow(self.inv(a), -n)
        result = 1
        while n:
            if n & 1:
                result = self.mul(result, a)
            a, n = self.mul(a, a), n >> 1
        return result

    def sigma(self, a):
        return self.elt(a % self.p, -(a // self.p))

    def sqrt(self, a):
        """Square root in K, or None if a is not a square in K."""
        p = self.p
        ar, ai = a % p, a // p
        if ai == 0:
            r = pow(ar, (p + 1) // 4, p)
            if r * r % p == ar:
                return r
            return self.elt(0, pow(-ar % p, (p + 1) // 4, p))
        norm = (ar * ar + ai * ai) % p
        r = pow(norm, (p + 1) // 4, p)
        if r * r % p != norm:
            return None
        for signed_r in (r, -r):
            u2 = (ar + signed_r) * pow(2, -1, p) % p
            u = pow(u2, (p + 1) // 4, p)
            if u and u * u % p == u2:
                root = self.elt(u, ai * pow(2 * u, -1, p))
                assert self.mul(root, root) == a
                return root
        raise ArithmeticError("square root invariant failed")

    def describe(self, a):
        return [str(a % self.p), str(a // self.p)]


class Fp4:
    """Only used to split the three quadratic factors of a Richelot neighbor."""

    def __init__(self, base):
        self.base, self.p = base, base.p
        self.q, self.size = base.size, base.size**2
        self.delta = base.p
        while base.sqrt(self.delta) is not None:
            self.delta = base.add(self.delta, 1)
        self.sigma_v_scale = base.pow(self.delta, (base.p - 1) // 2)

    def add(self, a, b):
        f, q = self.base, self.q
        return f.add(a % q, b % q) + q * f.add(a // q, b // q)

    def neg(self, a):
        return self.base.neg(a % self.q) + self.q * self.base.neg(a // self.q)

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        f, q = self.base, self.q
        ar, ai, br, bi = a % q, a // q, b % q, b // q
        real = f.add(f.mul(ar, br), f.mul(self.delta, f.mul(ai, bi)))
        imag = f.add(f.mul(ar, bi), f.mul(ai, br))
        return real + q * imag

    def inv(self, a):
        if not a:
            raise ZeroDivisionError
        f, q = self.base, self.q
        ar, ai = a % q, a // q
        scale = f.inv(f.sub(f.mul(ar, ar), f.mul(self.delta, f.mul(ai, ai))))
        return f.mul(ar, scale) + q * f.neg(f.mul(ai, scale))

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    pow = Fp2.pow

    def sigma(self, a):
        f, q = self.base, self.q
        return f.sigma(a % q) + q * f.mul(f.sigma(a // q), self.sigma_v_scale)

    def sqrt_base(self, a):
        root = self.base.sqrt(a)
        if root is not None:
            return root
        root = self.base.sqrt(self.base.div(a, self.delta))
        assert root is not None
        return self.q * root


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and not poly[-1]:
        poly.pop()
    return poly


def pmul(f, a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = f.add(out[i + j], f.mul(ai, bj))
    return trim(out)


def peval(f, poly, x):
    result = 0
    for c in reversed(poly):
        result = f.add(f.mul(result, x), c)
    return result


def bracket(f, a, b):
    da = [f.mul(i % f.p, a[i]) for i in range(1, len(a))] or [0]
    db = [f.mul(i % f.p, b[i]) for i in range(1, len(b))] or [0]
    u, v = pmul(f, da, b), pmul(f, a, db)
    return trim([f.sub(u[i] if i < len(u) else 0,
                       v[i] if i < len(v) else 0)
                 for i in range(max(len(u), len(v)))])


def determinant(f, rows):
    a, b, c = [list(row) + [0] * (3 - len(row)) for row in rows]
    return f.add(f.sub(f.mul(a[0], f.sub(f.mul(b[1], c[2]), f.mul(b[2], c[1]))),
                       f.mul(a[1], f.sub(f.mul(b[0], c[2]), f.mul(b[2], c[0])))),
                 f.mul(a[2], f.sub(f.mul(b[0], c[1]), f.mul(b[1], c[0]))))


def pairings(items):
    if not items:
        yield ()
        return
    a = items[0]
    for j in range(1, len(items)):
        for rest in pairings(items[1:j] + items[j + 1:]):
            yield ((a, items[j]),) + rest


def pair_poly(f, a, b):
    # Infinity is the missing sixth branch point of a quintic model.
    if a is None:
        return [f.neg(b), 1]
    if b is None:
        return [f.neg(a), 1]
    return [f.mul(a, b), f.neg(f.add(a, b)), 1]


def richelot(f, branches, partition):
    """Return determinant, three target factors, and the correctly twisted RHS.

    Assumes the source RHS is the monic product of its finite branch factors.
    Ascending coefficients: Delta*y^2 = [A,B][A,C][B,C].
    The scalar does not affect the branch-locus test, but matters for point counts.
    """
    a, b, c = [pair_poly(f, branches[i], branches[j]) for i, j in partition]
    delta = determinant(f, (a, b, c))
    if not delta:
        return 0, None, None
    factors = [bracket(f, a, b), bracket(f, a, c), bracket(f, b, c)]
    rhs = pmul(f, pmul(f, factors[0], factors[1]), factors[2])
    return delta, factors, [f.div(x, delta) for x in rhs]


def quadratic_roots(extension, poly):
    """Projective zeros of a binary quadratic; coefficients must lie in K."""
    f = extension.base
    poly = trim(poly)
    if len(poly) == 2:
        return [f.div(f.neg(poly[0]), poly[1]), None]
    if len(poly) != 3:
        raise ValueError("degenerate quadratic")
    c, b, a = poly
    disc = f.sub(f.mul(b, b), f.mul(4 % f.p, f.mul(a, c)))
    if not disc:
        raise ValueError("repeated branch point")
    root = extension.sqrt_base(disc)
    scale = f.inv(f.mul(2, a))
    return [extension.mul(extension.add(f.neg(b), root), scale),
            extension.mul(extension.sub(f.neg(b), root), scale)]


def point(x):
    return (1, 0) if x is None else (x, 1)


def cross(f, a, b):
    return f.sub(f.mul(a[0], b[1]), f.mul(b[0], a[1]))


def triple_matrix(f, triple):
    """Unique projective map taking triple to infinity, zero, one."""
    p0, p1, p2 = triple
    u, v = cross(f, p2, p0), cross(f, p2, p1)
    return (f.mul(u, p1[1]), f.neg(f.mul(u, p1[0])),
            f.mul(v, p0[1]), f.neg(f.mul(v, p0[0])))


def mmul(f, a, b):
    return (f.add(f.mul(a[0], b[0]), f.mul(a[1], b[2])),
            f.add(f.mul(a[0], b[1]), f.mul(a[1], b[3])),
            f.add(f.mul(a[2], b[0]), f.mul(a[3], b[2])),
            f.add(f.mul(a[2], b[1]), f.mul(a[3], b[3])))


def madj(f, a):
    return (a[3], f.neg(a[1]), f.neg(a[2]), a[0])


def mnormalize(f, a):
    scale = f.inv(next(x for x in a if x))
    return tuple(f.mul(x, scale) for x in a)


def mapply(f, a, p):
    return (f.add(f.mul(a[0], p[0]), f.mul(a[1], p[1])),
            f.add(f.mul(a[2], p[0]), f.mul(a[3], p[1])))


def branch_screen(extension, branches):
    """Necessary coordinate-invariant test for y^2=(x-alpha)*h(x), h in Fp[x].

    Such a model has a K-rational exceptional branch alpha, and the other
    five branches (including infinity) descend under a conjugate of p-Frobenius.
    Enumerate all possible maps using images of three distinct points. Positive
    results concern only branch data; the curve twist and an attack map remain.
    """
    f, k = extension, extension.base
    if len(branches) != 6 or len(set(branches)) != 6:
        raise ValueError("need six distinct geometric branch points")
    if any(x is not None and not 0 <= x < f.size for x in branches):
        raise ValueError("point outside splitting field representation")
    if {None if x is None else f.sigma(f.sigma(x)) for x in branches} != set(branches):
        raise ValueError("branch locus must be defined over K")
    witnesses, tried = [], 0
    for omitted, alpha in enumerate(branches):
        if alpha is not None and alpha >= k.size:
            continue  # the exceptional root must be K-rational
        t = [point(x) for i, x in enumerate(branches) if i != omitted]
        sigma_t = [tuple(f.sigma(c) for c in p) for p in t]
        source = triple_matrix(f, sigma_t[:3])
        seen = set()
        for targets in itertools.permutations(t, 3):
            tried += 1
            target = triple_matrix(f, targets)
            matrix = mnormalize(f, mmul(f, madj(f, target), source))
            if matrix in seen:
                continue
            seen.add(matrix)
            if any(c >= k.size for c in matrix):
                continue  # M must be defined over K, not just Fp^4
            cocycle = mmul(f, matrix, tuple(f.sigma(c) for c in matrix))
            if cocycle[1] or cocycle[2] or cocycle[0] != cocycle[3]:
                continue
            if all(any(cross(f, mapply(f, matrix, p), q) == 0 for q in t)
                   for p in sigma_t):
                witnesses.append({"omitted_branch_index": omitted,
                                  "matrix": [k.describe(c) for c in matrix]})
    return {"branch_condition_candidates": len(witnesses),
            "candidate_maps_tested": tried, "witnesses": witnesses,
            "status": "branch_candidate_needs_twist_and_map_checks" if witnesses
                      else "excluded_for_this_model_class"}


def audit_rosenhain(p, coefficients):
    k = Fp2(p)
    if len(coefficients) != 3:
        raise ValueError("need three Rosenhain parameters")
    for c in coefficients:
        validate_coefficient(p, c)
    f = Fp4(k)
    branches = [None, 0, 1] + [k.elt(*c) for c in coefficients]
    source = branch_screen(f, branches)
    neighbors = []
    for partition in pairings(tuple(range(6))):
        delta, factors, rhs = richelot(k, branches, partition)
        row = {"partition": partition, "determinant": k.describe(delta)}
        if not delta:
            row["status"] = "product_quotient_not_screened_as_genus2"
        else:
            roots = [root for factor in factors for root in quadratic_roots(f, factor)]
            row.update({"rhs_coefficients_ascending": [k.describe(x) for x in rhs],
                        "quadratic_factors_ascending": [[k.describe(x) for x in h]
                                                        for h in factors],
                        "K_rational_branch_points": sum(x is None or x < k.size
                                                        for x in roots),
                        "screen": branch_screen(f, roots),
                        "status": "screened"})
        neighbors.append(row)
    return {"p": str(p), "field": "Fp[i]/(i^2+1)",
            "coefficient_encoding": "[real, imaginary], decimal strings",
            "rosenhain_parameters": [[str(x) for x in c] for c in coefficients],
            "splitting_field_delta": k.describe(f.delta),
            "source": source, "neighbors": neighbors,
            "partitions_tested": len(neighbors),
            "nonsingular_genus2_quotients": sum(bool(x["determinant"] != ["0", "0"])
                                                for x in neighbors)}


def audit_bk(p, a):
    k = Fp2(p)
    validate_coefficient(p, a)
    a = k.elt(*a)
    if a == 0 or p == 5:
        raise ValueError("singular BK curve")
    if (k.size - 1) % 5:
        raise ValueError("test assumes fifth roots of unity in K")
    residue = k.pow(k.neg(a), (k.size - 1) // 5)
    irreducible = residue != 1 and a != 0
    return {"p": str(p), "a": k.describe(a),
            "minus_a_fifth_power_residue": k.describe(residue),
            "quintic_irreducible": irreducible,
            "branch_frobenius_cycle_type": [1, 5] if irreducible else [1] * 6,
            "rational_richelot_kernels": 0 if irreducible else 15,
            "reason": "A 5-cycle cannot preserve a perfect matching of six branches."
                      if irreducible else "All five finite roots are K-rational."}


def validate_coefficient(p, c):
    if not isinstance(c, (list, tuple)) or len(c) != 2:
        raise ValueError("coefficient must be [real, imaginary]")
    if any(type(x) is not int or not 0 <= x < p for x in c):
        raise ValueError("coefficient coordinates must be integers in [0,p)")


# Appendix A.2 and A.3 of the March 2013 Microsoft preprint (8glv.pdf).
# Each parameter is [real, imaginary], opposite to the PDF's printed ai+b order.
ROSENHAIN = [
    ("generic_montgomery", (2**31 - 307656) * 2**32 - 1,
     [(3027622460227942962, 2535046569008734284),
      (3306366075534497189, 6168887938445916986),
      (556271600778613627, 1342856250081051940)]),
    ("generic_nist_friendly", 2**64 - 189,
     [(18323491867119280217, 2952413239496295680),
      (17939661102861094359, 17185170913358942570),
      (11935307278311786128, 2480178925665143594)]),
    ("generic_mersenne", 2**61 - 1,
     [(895962277282143337, 819188565154922720),
      (1949052274256106659, 1814120660984670986),
      (1066444538207623350, 366557254231698774)]),
]
BK = [
    ("bk_montgomery", (2**31 - 201) * 2**32 - 1,
     (5793761506378637713, 4923475266616025179)),
    ("bk_nist_friendly", 2**64 - 2285,
     (9597225893510653377, 17275736705428938862)),
    ("bk_mersenne", 2**61 - 1, (1917550166465917276, 389685887377831749)),
]


def run_audit():
    return {
        "source": "https://www.microsoft.com/en-us/research/wp-content/uploads/2013/03/8glv.pdf",
        "source_pdf_sha256": "fc70bbcd657cd7a0f0c16671827f5eae3c19fd9f74adda49ea94630768071fe2",
        "parameter_sections": "Appendix A.2 (BK) and Appendix A.3 (Rosenhain/Kummer)",
        "scope": run_scope(),
        "bk": [{"name": name, **audit_bk(p, a)} for name, p, a in BK],
        "generic_and_kummer": [{"name": name, **audit_rosenhain(p, c)}
                               for name, p, c in ROSENHAIN],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        help="custom JSON: family, p, and parameters (Rosenhain) or a (BK)")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        if args.input:
            data = json.loads(args.input.read_text())
            if data['family'] == 'rosenhain':
                result = audit_rosenhain(data['p'], data['parameters'])
            elif data['family'] == 'bk':
                result = audit_bk(data['p'], data['a'])
            else:
                raise ValueError("family must be rosenhain or bk")
            result = {"scope": run_scope(), "input": data, "audit": result}
        else:
            result = run_audit()
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.error(str(error))
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


def run_scope():
    """Shared scope without rerunning the built-in audit."""
    return {"genus2_isogeny_depth": 1, "isogeny_type": "K-rational (2,2), Richelot",
            "branch_test_only": True, "DLP_solved": False,
            "security_certificate": False, "descent_map_constructed": False,
            "higher_degree_isogenies_tested": False,
            "novel_vulnerability_established": False,
            "target_model": "y^2=epsilon*(x-alpha)*h(x), alpha and epsilon in Fp2, h in Fp[x]"}


if __name__ == "__main__":
    main()
