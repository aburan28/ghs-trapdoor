"""Exact toy example: Richelot neighbor -> weak genus-2 model -> genus-4 cover.

No cryptographic-size DLP is solved. This constructs the cover morphism and
checks the complete Frobenius polynomial by exhaustive point counting.
Run: python -m research.genus2_cover_example --output research/results/toy-cover.json
"""

import argparse
import json
from pathlib import Path

from .genus2_audit import (Fp2, Fp4, branch_screen, pair_poly, pmul,
                           peval, quadratic_roots, richelot)


class Fp:
    def __init__(self, p):
        self.p = self.size = p

    def add(self, a, b):
        return (a + b) % self.p

    def neg(self, a):
        return -a % self.p

    def mul(self, a, b):
        return a * b % self.p


class F343:
    """F7[z]/(z^3-3); 3 is not a cube in F7, so the cubic is irreducible."""
    p, size = 7, 343

    @staticmethod
    def digits(a):
        return [a % 7, a // 7 % 7, a // 49]

    def add(self, a, b):
        return sum((x + y) % 7 * 7**i
                   for i, (x, y) in enumerate(zip(self.digits(a), self.digits(b))))

    def mul(self, a, b):
        x, y = self.digits(a), self.digits(b)
        c = [sum(x[j] * y[i-j] for j in range(3) if 0 <= i-j < 3)
             for i in range(5)]
        return (c[0] + 3*c[3]) % 7 + 7*((c[1] + 3*c[4]) % 7) + 49*(c[2] % 7)


def curve_count(field, rhs):
    squares = {field.mul(y, y) for y in range(field.size)}
    count = 0
    for x in range(field.size):
        value = peval(field, rhs, x)
        count += 1 if value == 0 else 2 if value in squares else 0
    degree = len(rhs) - 1
    return count + (1 if degree % 2 else 2 if rhs[-1] in squares else 0)


def poly_power(field, poly, n):
    out = [1]
    for _ in range(n):
        out = pmul(field, out, poly)
    return out


def squarefree_prime(poly, p):
    a, b = list(poly), [(i*poly[i]) % p for i in range(1, len(poly))]
    while b:
        while len(a) >= len(b):
            shift = len(a)-len(b)
            scale = a[-1] * pow(b[-1], -1, p) % p
            for i, coefficient in enumerate(b):
                a[i+shift] = (a[i+shift]-scale*coefficient) % p
            while a and a[-1] == 0:
                a.pop()
        a, b = b, a
    return len(a) == 1


def build_example():
    k = Fp2(7)
    extension = Fp4(k)
    base = Fp(7)
    # C: y^2=(x-(1+i))*h(x), h=x(x-1)(x-2)(x-3), i^2=-1.
    h, rhs = [1], [1]
    branches = [None, 0, 1, 2, 3, k.elt(1, 1)]
    for x in (0, 1, 2, 3):
        h = pmul(base, h, [-x % 7, 1])
    for x in branches[1:]:
        rhs = pmul(k, rhs, pair_poly(k, None, x))
    partition = ((0, 1), (2, 4), (3, 5))
    delta, factors, neighbor = richelot(k, branches, partition)
    neighbor_branches = [x for factor in factors for x in quadratic_roots(extension, factor)]
    screen = branch_screen(extension, neighbor_branches)
    assert screen['branch_condition_candidates'] == 0

    # x=N/(2t), N=-t^2+2t+1. Homogenize h to degree 4.
    numerator, denominator = [1, 2, 6], [0, 2]
    homogeneous = [0] * 9
    for j, coefficient in enumerate(h):
        term = pmul(base, poly_power(base, numerator, j),
                    poly_power(base, denominator, 4-j))
        for i, value in enumerate(term):
            homogeneous[i] = (homogeneous[i] + coefficient*value) % 7
    # H: S^2=2t*G(N,2t), with y=(1-i*t)*S/(8t^3).
    cover = pmul(base, [0, 2], homogeneous)
    assert squarefree_prime(cover, 7)
    # Verify the rational-function identity coefficient by coefficient:
    # (1-i*t)^2*H(t) = 2*t*F(N(t),2*t), with F the degree-5 homogenization.
    source_homogeneous = [0] * 11
    for j, coefficient in enumerate(rhs):
        term = pmul(k, poly_power(k, numerator, j),
                    poly_power(k, denominator, 5-j))
        for i, value in enumerate(term):
            source_homogeneous[i] = k.add(source_homogeneous[i], k.mul(coefficient, value))
    left = pmul(k, poly_power(k, [1, k.neg(k.p)], 2), cover)
    right = pmul(k, [0, 2], source_homogeneous)
    assert left == right

    n1, n2 = curve_count(k, rhs), curve_count(extension, rhs)
    c1 = n1 - k.size - 1
    c2 = (c1*c1 + n2 - k.size**2 - 1) // 2
    source_frob = [k.size**2, k.size*c1, c2, c1, 1]
    neighbor_counts = [curve_count(k, neighbor), curve_count(extension, neighbor)]
    assert neighbor_counts == [n1, n2]

    fields = (base, k, F343(), extension)
    counts = [curve_count(field, cover) for field in fields]
    powers = [7**n + 1 - counts[n-1] for n in range(1, 5)]
    descending = [1]
    for n in range(1, 5):
        value = -sum(descending[n-j]*powers[j-1] for j in range(1, n+1))
        assert value % n == 0
        descending.append(value // n)
    cover_frob = [0] * 9
    for i, coefficient in enumerate(descending):
        cover_frob[8-i] = coefficient
        cover_frob[i] = coefficient * 7**(4-i)
    expected = [0] * 9
    for i, coefficient in enumerate(source_frob):
        expected[2*i] = coefficient
    assert cover_frob == expected
    assert sum(source_frob) == sum(cover_frob) == 2**7 * 17

    checked = 0
    for t in range(1, k.size):
        root = k.sqrt(peval(k, cover, t))
        if root is None:
            continue
        for s in {root, k.neg(root)}:
            x = k.div(peval(k, numerator, t), k.mul(2, t))
            y = k.div(k.mul(k.sub(1, k.mul(k.p, t)), s), k.pow(t, 3))
            # 8=1 in characteristic 7, so denominator 8*t^3 is t^3 here.
            assert k.mul(y, y) == peval(k, rhs, x)
            checked += 1

    return {
        "scope": {"toy_only": True, "novel_attack_claimed": False,
                  "cover_map_constructed": True, "Richelot_point_map_implemented": False,
                  "DLP_solved": False, "practical_speedup_claimed": False},
        "field": "K=F7[i]/(i^2+1)",
        "coefficient_encoding": "[real, imaginary], ascending powers",
        "weak_curve": {"equation": "y^2=(x-(1+i))*x*(x-1)*(x-2)*(x-3)",
                       "rhs": [k.describe(x) for x in rhs]},
        "richelot_neighbor": {"rhs": [k.describe(x) for x in neighbor],
                              "determinant": k.describe(delta),
                              "partition_of_weak_curve": partition,
                              "branch_screen": screen,
                              "connection": "The dual Richelot isogeny maps this neighbor to the weak curve."},
        "cover": {"base_field": "F7", "genus": 4,
                  "equation": "S^2=" + "+".join(f"{c}*t^{i}" for i,c in enumerate(cover) if c),
                  "rhs": cover,
                  "squarefree_verified": True,
                  "rational_function_identity_verified": True,
                  "map_to_weak_curve": {"x": "(-t^2+2t+1)/(2t)",
                                        "y": "(1-i*t)*S/(8*t^3)"},
                  "ordinary_map_points_checked_over_F49": checked,
                  "exceptional_map_charts_not_evaluated": ["t=0", "t=infinity"]},
        "point_counts": {"weak_curve_F49_F2401": [n1, n2],
                         "neighbor_F49_F2401": neighbor_counts,
                         "cover_F7_F49_F343_F2401": counts},
        "frobenius_polynomials_ascending": {"weak_and_neighbor_over_F49": source_frob,
                                            "cover_over_F7": cover_frob,
                                            "cover_equals_source_evaluated_at_T_squared": True},
        "jacobian_group_order": sum(source_frob),
        "jacobian_group_order_factorization": {"2": 7, "17": 1},
        "subgroup_note": "The degree-4 Richelot isogeny and its dual preserve the order-17 subgroup. No pointwise Jacobian transport was run.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    encoded = json.dumps(build_example(), indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
