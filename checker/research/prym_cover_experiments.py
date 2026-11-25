"""Explicit genus-3 covers with elliptic quotient and genus-2 Prym factor.

For a squarefree monic quartic f with f(0)!=0, set
  C: y^2=f(x^2), E: y^2=f(u), D: v^2=u*f(u).
The maps C->E, (x,y)->(x^2,y), and C->D, (x,y)->(x^2,xy), give
Jac(C) isogenous to E x Jac(D). Jac(D) maps to Prym(C/E) with
2-primary kernel; odd-order subgroups survive. This is a known cover
mechanism, not a new DLP attack. We verify equations and full zeta factors.
"""

import argparse
import itertools
import json
from pathlib import Path

from .genus2_cover_example import squarefree_prime
from .small_fields import (SmallField, curve_count, evaluate, factors,
                           frobenius_polynomial, integer_poly_mul, points)


def build_case(p, f, fields=None):
    if len(f) != 5 or f[-1] != 1 or not f[0] or not squarefree_prime(f, p):
        raise ValueError('requires a squarefree monic quartic with nonzero constant term')
    if any(type(c) is not int or not 0 <= c < p for c in f):
        raise ValueError('quartic coefficients must be in [0,p)')
    fields = fields or [SmallField(p, n) for n in (1, 2, 3)]
    c = [0]*9
    c[::2] = f
    d = [0]+list(f)
    if not squarefree_prime(c, p) or not squarefree_prime(d, p):
        raise ArithmeticError('cover or genus-2 quotient is singular')
    counts = {'C_genus3': [], 'E_genus1': [], 'D_genus2': []}
    for k in fields:
        nc, ne, nd = (curve_count(k, rhs) for rhs in (c, f, d))
        if nc != ne+nd-k.size-1:
            raise ArithmeticError('point count factorization failed')
        for name, value in zip(counts, (nc, ne, nd)):
            counts[name].append(value)
    pc = frobenius_polynomial(p, counts['C_genus3'])
    pe = frobenius_polynomial(p, counts['E_genus1'][:1])
    pd = frobenius_polynomial(p, counts['D_genus2'][:2])
    if pc != integer_poly_mul(pe, pd):
        raise ArithmeticError('full Frobenius polynomial factorization failed')
    checked = 0
    k = fields[0]
    for x, y in points(k, c):
        u, v = k.mul(x, x), k.mul(x, y)
        if k.mul(y, y) != evaluate(k, f, u) or k.mul(v, v) != evaluate(k, d, u):
            raise ArithmeticError('cover point map failed')
        checked += 1
    ne, nd, nc = sum(pe), sum(pd), sum(pc)
    if nc != ne*nd:
        raise ArithmeticError('Jacobian order factorization failed')
    prym_primes = [r for r, _ in factors(nd) if r != 2 and ne % r]
    return {'p': p, 'f_ascending': f, 'C_rhs_ascending': c, 'D_rhs_ascending': d,
            'maps': {'C_to_E': ['u=x^2', 'y=y'], 'C_to_D': ['u=x^2', 'v=x*y']},
            'point_counts_over_p_p2_p3': counts,
            'frobenius_polynomials_ascending': {'C': pc, 'E': pe, 'D': pd},
            'jacobian_orders': {'C': nc, 'E': ne, 'D': nd},
            'D_order_factorization': factors(nd),
            'odd_primes_in_D_order_absent_from_E_order': prym_primes,
            'rational_affine_C_points_map_checked': checked,
            'both_points_at_infinity_included_in_C_count': True,
            'symbolic_map_identities': ['y^2=f(x^2)', '(x*y)^2=x^2*f(x^2)'],
            'subgroup_preservation_reason': 'pi_D_pushforward o pi_D_pullback = [2]',
            'pointwise_Jacobian_divisor_arithmetic_implemented': False}


def run():
    cases = []
    for p in (5, 7, 11, 13):
        fields = [SmallField(p, n) for n in (1, 2, 3)]
        found = 0
        # Exact bounded family: first eight squarefree quartics in this order.
        for constant, linear, quadratic, cubic in itertools.product(range(p), repeat=4):
            f = [constant, linear, quadratic, cubic, 1]
            if not constant or not squarefree_prime(f, p):
                continue
            cases.append(build_case(p, f, fields))
            found += 1
            if found == 8:
                break
    return {'scope': {'toy_only': True, 'known_cover_mechanism': True,
                      'Prym_dimension': 2, 'Prym_isogenous_to_genus2_Jacobian': True,
                      'explicit_curve_maps_constructed': True,
                      'complete_Frobenius_polynomials_checked': True,
                      'subgroup_injectivity_for_odd_order_follows_from_degree2_adjunction': True,
                      'infinity_map_charts_explicitly_evaluated': False,
                      'Jacobian_DLP_solved': False, 'practical_speedup_established': False,
                      'novel_vulnerability_established': False},
            'search_bounds': {'base_primes': [5, 7, 11, 13], 'quartics_per_prime': 8,
                              'ordering': '(constant,linear,quadratic,cubic) lexicographic'},
            'case_count': len(cases),
            'cases_with_odd_Prym_prime_absent_from_elliptic_factor': sum(
                bool(c['odd_primes_in_D_order_absent_from_E_order']) for c in cases),
            'cases': cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result)
    else:
        print(result, end='')


if __name__ == '__main__':
    main()
