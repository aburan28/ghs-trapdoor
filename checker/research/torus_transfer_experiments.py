"""Positive/negative controls for descent to a norm-one torus.

The smooth locus of the singular nodal curve y^2=x^2(x+delta) is
isomorphic to ker(N:Fp2* -> Fp*) for nonsquare delta. We verify its full
group law and toy DLP transfer. The same rational function fails the group
law on the smooth elliptic perturbation y^2=x^3+delta*x^2+1.
"""

import argparse
import itertools
import json
from pathlib import Path

from .small_fields import SmallField, ec_add, ec_mul, factors, points


def experiment(p):
    k = SmallField(p)
    delta = next(x for x in range(2, p) if pow(x, (p-1)//2, p) == p-1)
    extension = SmallField(p, 2, [-delta % p, 0, 1])
    gamma = p  # gamma^2=delta
    group = [None]+list(points(k, [0, 0, delta, 1]))
    group.remove((0, 0))  # the node is not in the group

    def to_torus(point):
        if point is None:
            return 1
        x, y = point
        return extension.div(extension.sub(y, extension.mul(gamma, x)),
                             extension.add(y, extension.mul(gamma, x)))

    def from_torus(t):
        if t == 1:
            return None
        z = extension.div(extension.mul(gamma, extension.add(1, t)), extension.sub(1, t))
        if z >= p:
            raise ArithmeticError('inverse parameter did not descend to Fp')
        x = (z*z-delta) % p
        return x, z*x % p

    norm_one = {t for t in range(1, extension.size) if extension.pow(t, p+1) == 1}
    image = {to_torus(point) for point in group}
    if len(group) != p+1 or image != norm_one:
        raise ArithmeticError('torus bijection failed')
    if not all(from_torus(to_torus(point)) == point for point in group):
        raise ArithmeticError('inverse torus map failed')
    for left, right in itertools.product(group, repeat=2):
        if to_torus(ec_add(k, delta, 0, left, right)) != \
                extension.mul(to_torus(left), to_torus(right)):
            raise ArithmeticError('nodal chord law does not match torus multiplication')
    generator = next(point for point in group if all(
        ec_mul(k, delta, 0, point, (p+1)//r) is not None for r, _ in factors(p+1)))
    g = to_torus(generator)
    logs = {extension.pow(g, n): n for n in range(p+1)}
    for n in range(p+1):
        if logs[to_torus(ec_mul(k, delta, 0, generator, n))] != n:
            raise ArithmeticError('DLP transfer failed')
    secret = (p+1)//2+1
    target = ec_mul(k, delta, 0, generator, secret)

    # An ordinary elliptic curve with the same quadratic coefficient.
    smooth = [None]+list(points(k, [1, 0, delta, 1]))
    # For y²=x³+delta*x²+1, nonsingularity is 4*delta³+27 != 0.
    if (4*delta**3+27) % p == 0:
        raise ArithmeticError('negative control must be smooth')
    counterexample = None
    for left, right in itertools.product(smooth, repeat=2):
        total = ec_add(k, delta, 0, left, right)
        product_image = extension.mul(to_torus(left), to_torus(right))
        if to_torus(total) != product_image:
            counterexample = {'P': left, 'Q': right, 'P_plus_Q': total,
                              'map_P_plus_Q': to_torus(total),
                              'map_P_times_map_Q': product_image}
            break
    if counterexample is None:
        raise ArithmeticError('negative control did not refute the candidate map')
    return {'p': p, 'delta_nonsquare': delta,
            'quadratic_modulus_ascending': extension.modulus,
            'nodal_equation': 'y^2=x^2*(x+delta)',
            'torus_map': '(y-gamma*x)/(y+gamma*x); gamma^2=delta; O maps to 1',
            'inverse_map': 'z=gamma*(1+t)/(1-t); x=z^2-delta; y=z*x; t=1 maps to O',
            'smooth_locus_order': len(group), 'torus_order': len(norm_one),
            'bijection_and_inverse_verified': True,
            'group_law_pairs_checked': len(group)**2,
            'all_scalar_logs_checked': p+1,
            'toy_DLP': {'generator': generator, 'target': target, 'generator_torus': g,
                        'target_torus': to_torus(target), 'secret': secret,
                        'recovered_by_torus_enumeration': logs[to_torus(target)]},
            'smooth_elliptic_negative_control': {
                'equation': 'y^2=x^3+delta*x^2+1', 'group_order': len(smooth),
                'nonsingularity_verified': True,
                'same_candidate_map_is_homomorphism': False,
                'explicit_counterexample': counterexample}}


def run():
    return {'scope': {'toy_only': True, 'known_singular_curve_to_torus_isomorphism': True,
                      'DLP_transfer_demonstrated_on_singular_curves_only': True,
                      'smooth_elliptic_to_torus_reduction_established': False,
                      'novel_vulnerability_established': False},
            'experiments': [experiment(p) for p in (5, 7, 11, 13)]}


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
