"""Actual small-curve census for rational 2/3-torsion quotient kernels.

Exhaust all short Weierstrass equations over F125; deterministically sample
512, 64 and 64 nonsingular equations over F343, F1331 and F2197. Equations
are not deduplicated by isomorphism class. No quotient target is reconstructed.
"""

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

from .small_fields import SmallField, ec_mul, factors, points
from .weil_quotient_kernels import screen


def curve_data(k, a, b, squares=None, powers=None):
    """Count E and its rational division roots; determine Frobenius conjugacy classes."""
    p, size = k.p, k.size
    add, mul = k.add, k.mul
    a2, a3, b2 = mul(a, a), mul(mul(a, a), a), mul(b, b)
    discriminant = add(mul(4 % p, a3), mul(27 % p, b2))
    if discriminant == 0:
        return None
    if squares is None:
        squares = {mul(y, y) for y in range(size)}
    if powers is None:
        powers = [(mul(x, x), mul(mul(x, x), x), mul(mul(x, x), mul(x, x)))
                  for x in range(size)]
    negative_a2 = k.neg(a2)
    count, two_roots, three_roots, rational_three_points = 1, [], [], 0
    for x, (x2, x3, x4) in enumerate(powers):
        value = add(add(x3, mul(a, x)), b)
        count += 1 if value == 0 else 2 if value in squares else 0
        if value == 0:
            two_roots.append(x)
        psi3 = add(add(mul(3, x4), mul(6 % p, mul(a, x2))),
                   add(mul(12 % p, mul(b, x)), negative_a2))
        if psi3 == 0:
            three_roots.append(x)
            if value in squares:
                rational_three_points += 2
    trace = size+1-count
    if trace*trace > 4*size:
        raise ArithmeticError('Hasse bound failed')
    if len(two_roots) == 3:
        m2 = ((1, 0), (0, 1))
    elif len(two_roots) == 1:
        m2 = ((0, 1), (1, 0))
    elif len(two_roots) == 0:
        m2 = ((0, 1), (1, 1))
    else:
        raise ArithmeticError('unexpected 2-division root count')
    t, determinant = trace % 3, size % 3
    repeated = (t*t-4*determinant) % 3 == 0
    if repeated and len(three_roots) == 4:
        eigenvalue = 2*t % 3
        m3 = ((eigenvalue, 0), (0, eigenvalue))
    else:
        if repeated and len(three_roots) != 1:
            raise ArithmeticError('repeated-root Frobenius classification failed')
        m3 = ((0, 1), (-determinant % 3, t))
    j = k.div(mul(1728 % p, mul(4 % p, a3)), discriminant)
    return {'a': a, 'b': b, 'group_order': count, 'trace': trace,
            'two_division_roots': two_roots, 'three_division_roots': three_roots,
            'rational_three_torsion_points_including_identity': 1+rational_three_points,
            'frobenius_mod_2_representative': m2, 'frobenius_mod_3_representative': m3,
            'j_in_prime_field': k.pow(j, p) == j, 'j': j}


@lru_cache(maxsize=None)
def counts_for_module(q, ell, matrix):
    report = screen(q, ell, matrix, witness_limit=0)
    return (report['frobenius_stable_lagrangians'],
            report['stable_kernels_passing_no_coordinate_factor_test'])


def sample_pairs(k, count):
    """SHA256 counter stream, no replacement; only nonsingular pairs count toward bounds."""
    seen, counter = set(), 0
    while len(seen) < k.size**2:
        digest = sha256(f'ghs-cubic-census-v1:{k.p}:{counter}'.encode()).digest()
        counter += 1
        pair = (int.from_bytes(digest[:16], 'big') % k.size,
                int.from_bytes(digest[16:], 'big') % k.size)
        if pair not in seen:
            seen.add(pair)
            yield pair
            count -= 1
            if count == 0:
                return


def subgroup_witness(k, row):
    """Find an actual nonzero point of the largest prime order coprime to 6, if possible."""
    n, a, b = row['group_order'], row['a'], row['b']
    eligible = [r for r, _ in factors(n) if r > 3]
    for r in reversed(eligible):
        for point in points(k, [b, a, 0, 1]):
            candidate = ec_mul(k, 0, a, point, n//r)
            if candidate is not None and ec_mul(k, 0, a, candidate, r) is None:
                return {'order': r, 'point': candidate,
                        'order_verified_by_scalar_multiplication': True,
                        'coprime_to_2_and_3': True}
    return None


def census(p, sample_size=None):
    k = SmallField(p, 3)
    if k.size <= 343:
        k.tables()
    squares = {k.mul(y, y) for y in range(k.size)}
    powers = [(k.pow(x, 2), k.pow(x, 3), k.pow(x, 4)) for x in range(k.size)]
    pairs = ((a, b) for a in range(k.size) for b in range(k.size)) if sample_size is None else \
            sample_pairs(k, min(k.size**2, 4*sample_size))
    totals, classes, representatives, digest = Counter(), Counter(), {}, sha256()
    prime_order_candidates = []
    for a, b in pairs:
        totals['coefficient_pairs_examined'] += 1
        row = curve_data(k, a, b, squares, powers)
        if row is None:
            totals['singular_equations_skipped'] += 1
            continue
        n = row['group_order']
        m2, m3 = row['frobenius_mod_2_representative'], row['frobenius_mod_3_representative']
        stable2, mixed2 = counts_for_module(p, 2, m2)
        stable3, mixed3 = counts_for_module(p, 3, m3)
        two_blocked_three_candidate = stable2 == 0 and mixed3 > 0
        n_factors = factors(n)
        prime_order = n_factors == [[n, 1]]
        totals['nonsingular_equations'] += 1
        totals['j_outside_prime_field'] += not row['j_in_prime_field']
        totals['with_mixed_2_kernel_candidates'] += mixed2 > 0
        totals['with_mixed_3_kernel_candidates'] += mixed3 > 0
        totals['zero_rational_2_kernels_but_mixed_3_candidates'] += two_blocked_three_candidate
        totals['prime_order_equations'] += prime_order
        totals['prime_order_zero_2_but_mixed_3_candidates'] += prime_order and two_blocked_three_candidate
        totals['prime_order_zero_2_mixed_3_and_j_outside_base'] += \
            prime_order and two_blocked_three_candidate and not row['j_in_prime_field']
        key = (m2, m3)
        classes[key] += 1
        record_key = f'module:{m2}:{m3}'
        prime_candidate = two_blocked_three_candidate and prime_order and not row['j_in_prime_field']
        if record_key not in representatives or prime_candidate:
            record = {
                **row, 'group_order_factorization': n_factors,
                'rational_2_kernels': stable2, 'mixed_2_kernel_candidates': mixed2,
                'rational_3_kernels': stable3, 'mixed_3_kernel_candidates': mixed3,
                'subgroup_witness': subgroup_witness(k, row)}
            representatives.setdefault(record_key, record)
            if prime_candidate:
                prime_order_candidates.append(record)
        digest.update((json.dumps([a, b, n, row['two_division_roots'],
                                  row['three_division_roots'], m2, m3],
                                 separators=(',', ':'))+'\n').encode())
        if sample_size is not None and totals['nonsingular_equations'] == sample_size:
            break
    if sample_size is not None and totals['nonsingular_equations'] != sample_size:
        raise ArithmeticError('sample stream exhausted before reaching bound')
    return {'p': p, 'extension_degree': 3, 'field_size': k.size,
            'modulus_ascending': k.modulus,
            'equation': 'y^2=x^3+a*x+b',
            'sampling': 'exhaustive coefficient pairs' if sample_size is None else
                        'SHA256 ghs-cubic-census-v1 counter stream; without replacement',
            'requested_nonsingular_sample_size': sample_size,
            'totals': dict(sorted(totals.items())), 'ordered_curve_data_sha256': digest.hexdigest(),
            'module_classes': [
                {'frobenius_mod_2_representative': m2, 'frobenius_mod_3_representative': m3,
                 'equations': count, 'two_kernel_counts': counts_for_module(p, 2, m2),
                 'three_kernel_counts': counts_for_module(p, 3, m3)}
                for (m2, m3), count in sorted(classes.items())],
            'prime_order_extension_j_candidates': prime_order_candidates,
            'representatives': representatives}


def run():
    return {'scope': {'toy_parameters_only': True, 'equations_not_isomorphism_classes': True,
                      'point_counts_exhaustive_for_each_equation': True,
                      'torsion_Frobenius_conjugacy_classes_verified': True,
                      'torsion_basis_coordinates_constructed': False,
                      'quotient_targets_constructed': False, 'DLP_solved': False,
                      'novel_vulnerability_established': False},
            'fields': [census(5), census(7, 512), census(11, 64), census(13, 64)]}


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
