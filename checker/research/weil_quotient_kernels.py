"""Screen rational polarized quotients of Res_{Fq^3/Fq}(E) on 2/3-torsion.

This enumerates kernels, not curve equations or DLP attacks. An explicit
Frobenius matrix is required; trace/determinant alone can be ambiguous.
Run: python -m research.weil_quotient_kernels --output /tmp/quotients.json
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import itertools
import json
from math import gcd
from pathlib import Path


def row_times(v, matrix, ell):
    return tuple(sum(x * matrix[i][j] for i, x in enumerate(v)) % ell
                 for j in range(len(matrix[0])))


def rank(rows, ell):
    if not rows:
        return 0
    a = [list(row) for row in rows]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j] % ell), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = pow(a[r][j] % ell, -1, ell)
        a[r] = [x * scale % ell for x in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][j] % ell
            a[i] = [(x - scale*y) % ell for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def rref_planes(ell):
    """Every 3-plane in F_ell^6 exactly once, in reduced row-echelon form."""
    for pivots in itertools.combinations(range(6), 3):
        free = [(i, j) for i in range(3) for j in range(pivots[i] + 1, 6)
                if j not in pivots]
        for values in itertools.product(range(ell), repeat=len(free)):
            rows = [[int(j == pivots[i]) for j in range(6)] for i in range(3)]
            for (i, j), value in zip(free, values):
                rows[i][j] = value
            yield tuple(tuple(row) for row in rows)


def pairing(u, v, ell, q):
    return sum(pow(q, i, ell) * (u[2*i]*v[2*i+1] - u[2*i+1]*v[2*i])
               for i in range(3)) % ell


def isotropic(rows, ell, q):
    return all(pairing(u, v, ell, q) == 0 for u, v in itertools.combinations(rows, 2))


def axis_intersections(rows, ell):
    # dim(L intersect V_i) = dim(L) - rank(projection to the other two factors).
    return tuple(3-rank([tuple(x for j, x in enumerate(row) if j//2 != i)
                         for row in rows], ell) for i in range(3))


@lru_cache(maxsize=4)
def lagrangians(ell, q_residue):
    if ell not in (2, 3) or q_residue % ell == 0:
        raise ValueError('only ell=2 or 3, coprime to q, is supported')
    result = tuple((rows, axis_intersections(rows, ell))
                   for rows in rref_planes(ell) if isotropic(rows, ell, q_residue))
    expected = (ell+1)*(ell**2+1)*(ell**3+1)
    if len(result) != expected:
        raise ArithmeticError('Lagrangian enumeration count failed')
    return result


def coordinates(v, rows, ell):
    pivots = [next(j for j, x in enumerate(row) if x) for row in rows]
    coefficients = tuple(v[j] for j in pivots)
    remainder = tuple((x-sum(coefficients[i]*rows[i][j] for i in range(3))) % ell
                      for j, x in enumerate(v))
    return coefficients if not any(remainder) else None


def frobenius_matrix(q, ell, matrix):
    """Row convention; bases on factors 1/2 are q/q^2-powers of factor 0."""
    if type(q) is not int or q < 2 or ell not in (2, 3) or gcd(q, ell) != 1:
        raise ValueError('requires q>=2 and ell=2 or 3 coprime to q')
    if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
        raise ValueError('need an explicit 2x2 q^3-Frobenius matrix')
    if any(type(x) is not int or not 0 <= x < ell for row in matrix for x in row):
        raise ValueError('matrix entries must be integers in [0,ell)')
    if (matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]) % ell != pow(q, 3, ell):
        raise ValueError('matrix determinant must equal q^3 modulo ell')
    f = [[0]*6 for _ in range(6)]
    f[0][2] = f[1][3] = f[2][4] = f[3][5] = 1
    f[4][0:2], f[5][0:2] = list(matrix[0]), list(matrix[1])
    return tuple(tuple(row) for row in f)


def screen(q, ell, matrix, subgroup_order=None, witness_limit=3):
    f = frobenius_matrix(q, ell, matrix)
    if subgroup_order is not None and (type(subgroup_order) is not int or subgroup_order < 2):
        raise ValueError('claimed subgroup order must be an integer >=2')
    if type(witness_limit) is not int or not 0 <= witness_limit <= 100:
        raise ValueError('witness_limit must be between 0 and 100')
    stable = mixed = 0
    witnesses = []
    for rows, intersections in lagrangians(ell, q % ell):
        action = [coordinates(row_times(row, f, ell), rows, ell) for row in rows]
        if any(v is None for v in action):
            continue
        stable += 1
        if any(intersections):
            continue
        mixed += 1
        if len(witnesses) < witness_limit:
            witnesses.append({'basis_rows': rows, 'frobenius_on_kernel_rows': action,
                              'coordinate_factor_intersection_dimensions': intersections})
    return {
        'q': q, 'ell': ell, 'extension_degree': 3, 'dimension': 3,
        'q3_frobenius_on_E_torsion_rows': matrix,
        'q_frobenius_on_Weil_restriction_rows': f,
        'pairing_block_weights': [pow(q, i, ell) for i in range(3)],
        'geometric_lagrangians_enumerated': len(lagrangians(ell, q % ell)),
        'frobenius_stable_lagrangians': stable,
        'stable_kernels_passing_no_coordinate_factor_test': mixed,
        'witnesses': witnesses,
        'quotient_isogeny_degree': ell**3,
        'claimed_subgroup_order': subgroup_order,
        'subgroup_injective_if_claimed_order_is_valid':
            gcd(subgroup_order, ell) == 1 if subgroup_order is not None else None,
        'status': 'kernel_candidates_require_target_and_map_construction' if mixed else
                  'no_kernel_passes_these_conditions_for_this_torsion_module',
        'scope': {'explicit_frobenius_matrix_required': True,
                  'q_prime_power_assumed': True,
                  'curve_realization_of_input_module_verified': False,
                  'geometric_indecomposability_verified': False,
                  'jacobian_over_base_field_verified': False,
                  'quotient_equations_constructed': False,
                  'point_transport_implemented': False,
                  'DLP_solved': False, 'novel_vulnerability_established': False},
    }


class Cubic523:
    """Fp[w]/(w^3+2), used only for the published example's exact 3-torsion check."""
    p, size = 523, 523**3

    def elt(self, a, b=0, c=0):
        return a % self.p + self.p*(b % self.p) + self.p**2*(c % self.p)

    def digits(self, x):
        return (x % self.p, x//self.p % self.p, x//self.p**2)

    def add(self, a, b):
        return self.elt(*(x+y for x, y in zip(self.digits(a), self.digits(b))))

    def neg(self, a):
        return self.elt(*(-x for x in self.digits(a)))

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        x, y = self.digits(a), self.digits(b)
        c = [sum(x[j]*y[i-j] for j in range(3) if 0 <= i-j < 3) for i in range(5)]
        return self.elt(c[0]-2*c[3], c[1]-2*c[4], c[2])

    def pow(self, a, n):
        out = 1
        while n:
            if n & 1:
                out = self.mul(out, a)
            a, n = self.mul(a, a), n >> 1
        return out

    def inv(self, a):
        if not a:
            raise ZeroDivisionError
        return self.pow(a, self.size-2)


def poly_strip(poly):
    poly = list(poly)
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def poly_remainder(k, a, b):
    a, b = poly_strip(a), poly_strip(b)
    while len(a) >= len(b) and a:
        shift = len(a)-len(b)
        scale = k.mul(a[-1], k.inv(b[-1]))
        for i, x in enumerate(b):
            a[i+shift] = k.sub(a[i+shift], k.mul(scale, x))
        a = poly_strip(a)
    return a


def poly_pow_mod(k, a, n, modulus):
    def mul(u, v):
        out = [0]*(len(u)+len(v)-1)
        for i, x in enumerate(u):
            for j, y in enumerate(v):
                out[i+j] = k.add(out[i+j], k.mul(x, y))
        return poly_remainder(k, out, modulus)
    out = [1]
    while n:
        if n & 1:
            out = mul(out, a)
        a, n = mul(a, a), n >> 1
    return out


def poly_gcd_degree(k, a, b):
    while b:
        a, b = b, poly_remainder(k, a, b)
    return len(poly_strip(a))-1


def verify_published_torsion():
    """Verify a basis of E[3] from Tian §3.4 without trusting a scalar trace guess."""
    k = Cubic523()
    if pow(-2 % 523, (523-1)//3, 523) == 1:
        raise ArithmeticError('cubic must be irreducible')
    a, b = k.elt(140, 92, 82), k.elt(243, 339, 359)
    delta = k.elt(6, 1)  # theta^2=w+6
    if k.pow(delta, (k.size-1)//2) != k.neg(1):
        raise ArithmeticError('quadratic must be irreducible')
    points = [(k.elt(337, 370, 222), k.elt(400, 387, 122)),
              (k.elt(156, 484, 319), k.elt(108, 484, 298))]
    for x, y_coefficient in points:
        y2 = k.mul(k.mul(y_coefficient, y_coefficient), delta)
        rhs = k.add(k.add(k.pow(x, 3), k.mul(a, x)), b)
        psi3 = k.sub(k.add(k.add(k.mul(3, k.pow(x, 4)),
                                 k.mul(6, k.mul(a, k.pow(x, 2)))),
                           k.mul(12, k.mul(b, x))), k.mul(a, a))
        if y2 != rhs or psi3 != 0 or not y_coefficient:
            raise ArithmeticError('published 3-torsion point verification failed')
    if points[0][0] == points[1][0]:
        raise ArithmeticError('need independent 3-torsion points')
    two_division = [b, a, 0, 1]
    remainder = poly_pow_mod(k, [0, 1], k.size, two_division)
    remainder += [0]*(2-len(remainder))
    remainder[1] = k.sub(remainder[1], 1)
    root_gcd_degree = poly_gcd_degree(k, two_division, poly_strip(remainder))
    if root_gcd_degree != 0:
        raise ArithmeticError('published 2-division polynomial should be irreducible')
    # Distinct x coordinates imply independence on E[3]. Both have x in K and
    # y=coefficient*theta, so q^3-Frobenius sends both to their negatives.
    return {'source': 'https://arxiv.org/html/2012.07173v1#S3.SS4',
            'q': 523, 'field': 'F523[w]/(w^3+2)',
            'theta_squared': [6, 1, 0],
            'basis_points_x_and_y_over_theta':
                [[k.digits(x), k.digits(y)] for x, y in points],
            'curve_membership_verified': True, 'three_division_polynomial_verified': True,
            'independence_verified_by_distinct_x': True,
            'q3_frobenius_matrix_verified': [[2, 0], [0, 2]],
            'two_division_cubic_irreducible_verified': True,
            'gcd_degree_with_X_to_field_order_minus_X': root_gcd_degree,
            'claimed_curve_group_order': 143069629,
            'curve_group_order_recounted': False}


def residue_survey():
    results = []
    for ell in (2, 3):
        for q_residue in range(1, ell):
            # Determinant of the E/Fq^3 Frobenius is q^3 mod ell.
            delta = pow(q_residue, 3, ell)
            for trace in range(ell):
                matrices = [('companion', [[0, 1], [-delta % ell, trace]])]
                matrices += [('scalar', [[lam, 0], [0, lam]]) for lam in range(1, ell)
                             if 2*lam % ell == trace and lam*lam % ell == delta]
                for kind, matrix in matrices:
                    row = screen(q_residue+ell, ell, matrix, witness_limit=1)
                    results.append({'matrix_type': kind, 'trace_residue': trace,
                                    'q_residue': q_residue, **row})
    return results


def run():
    verification = verify_published_torsion()
    report = screen(523, 3, verification['q3_frobenius_matrix_verified'],
                    subgroup_order=143069629)
    report['scope']['curve_realization_of_input_module_verified'] = True
    two_report = screen(523, 2, [[0, 1], [1, 1]], subgroup_order=143069629)
    two_report['scope']['curve_realization_of_input_module_verified'] = True
    return {'published_example_verification': verification,
            'published_example_kernel_screen': report,
            'same_curve_two_torsion_obstruction': two_report,
            'residue_class_survey_not_a_curve_vulnerability_survey': residue_survey()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        help='JSON: q, ell, frobenius_matrix, optional subgroup_order')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        if args.input:
            data = json.loads(args.input.read_text())
            report = screen(data['q'], data['ell'], data['frobenius_matrix'],
                            data.get('subgroup_order'))
        else:
            report = run()
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.error(str(error))
    encoded = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
