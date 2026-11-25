import itertools
import json
from pathlib import Path
import subprocess
import sys
import unittest

from research.cubic_curve_census import census, curve_data
from research.genus2_audit import Fp2
from research.genus2_cover_example import F343, squarefree_prime
from research.prym_cover_experiments import build_case, run as run_pryms
from research.small_fields import (SmallField, curve_count, ec_mul, evaluate,
                                   factors, points)
from research.torus_transfer_experiments import run as run_tori
from research.weil_quotient_kernels import rank


ROOT = Path(__file__).resolve().parents[1]


def saved(name):
    return json.loads((ROOT/'research/results'/f'{name}.json').read_text())


def rational_three_torsion_by_group_law(k, a, b):
    # Independent of division-polynomial roots or Frobenius classification.
    return 1+sum(ec_mul(k, 0, a, point, 3) is None for point in points(k, [b, a, 0, 1]))


class ExperimentBatchTests(unittest.TestCase):
    def test_small_fields_against_existing_independent_representations(self):
        k, old = SmallField(7, 2, [1, 0, 1]), Fp2(7)
        for a, b in itertools.product(range(49), repeat=2):
            self.assertEqual(k.add(a, b), old.add(a, b))
            self.assertEqual(k.mul(a, b), old.mul(a, b))
        k, old = SmallField(7, 3, [4, 0, 0, 1]), F343()
        for a in range(343):
            b = (a*23+41) % 343
            self.assertEqual(k.mul(a, b), old.mul(a, b))
            if a:
                self.assertEqual(k.mul(a, k.inv(a)), 1)
        before = [k.mul(a, 7) for a in range(343)]
        k.tables()
        self.assertEqual(before, [k.mul(a, 7) for a in range(343)])

    def test_small_field_validation_and_frobenius(self):
        for p, degree, modulus in ((9, 2, None), (3, 2, None), (7, 4, None),
                                    (7, 3, [1, 0, 0, 1]), (7, 2, [1, 1, 2])):
            with self.assertRaises(ValueError):
                SmallField(p, degree, modulus)
        for p in (5, 7):
            k = SmallField(p, 3)
            self.assertTrue(all(k.pow(a, k.size) == a for a in range(k.size)))
        with self.assertRaises(ValueError):
            SmallField(11, 3).tables()

    def test_derivative_degree_drop_and_infinity_counts(self):
        # x^5+x+1 has derivative 1 in characteristic 5; x^5+1 has derivative 0.
        self.assertTrue(squarefree_prime([1, 1, 0, 0, 0, 1], 5))
        self.assertFalse(squarefree_prime([1, 0, 0, 0, 0, 1], 5))
        k = SmallField(5)
        for f in ([1, 1, 0, 0, 1], [1, 1, 0, 0, 2], [1, 1, 0, 0, 0, 1]):
            affine = sum(1 for x in range(5) for y in range(5)
                         if y*y % 5 == evaluate(k, f, x))
            infinity = 1 if (len(f)-1) % 2 else 2 if f[-1] == 1 else 0
            self.assertEqual(curve_count(k, f), affine+infinity)

    def test_curve_census_representatives_recompute_exactly(self):
        for field in saved('cubic-curve-census')['fields']:
            k = SmallField(field['p'], 3, field['modulus_ascending'])
            if k.size <= 343:
                k.tables()
            for row in field['representatives'].values():
                actual = curve_data(k, row['a'], row['b'])
                for key in actual:
                    self.assertEqual(json.loads(json.dumps(actual[key])), row[key])
                witness = row['subgroup_witness']
                if witness:
                    r, point = witness['order'], tuple(witness['point'])
                    self.assertEqual(factors(r), [[r, 1]])
                    self.assertIsNotNone(point)
                    self.assertIsNone(ec_mul(k, 0, row['a'], point, r))

    def test_repeated_eigenvalues_verified_using_twist_group_law(self):
        field = saved('cubic-curve-census')['fields'][1]
        k = SmallField(7, 3, field['modulus_ascending'])
        k.tables()
        delta = next(x for x in range(1, k.size) if k.pow(x, (k.size-1)//2) != 1)
        observed_minus_identity_ranks = set()
        for row in field['representatives'].values():
            a, b = row['a'], row['b']
            m = row['frobenius_mod_3_representative']
            for sign in (1, -1):
                matrix = [[(m[i][j]-(sign if i == j else 0)) % 3 for j in range(2)]
                          for i in range(2)]
                predicted = 3**(2-rank(matrix, 3))
                aa, bb = (a, b) if sign == 1 else \
                         (k.mul(a, k.pow(delta, 2)), k.mul(b, k.pow(delta, 3)))
                actual = rational_three_torsion_by_group_law(k, aa, bb)
                self.assertEqual(actual, predicted)
                if sign == -1 and (m[0][0]+m[1][1]) % 3 == 1:
                    observed_minus_identity_ranks.add(actual)
        self.assertEqual(observed_minus_identity_ranks, {3, 9})

    def test_census_bounds_and_prime_order_candidates(self):
        report = saved('cubic-curve-census')
        self.assertFalse(report['scope']['quotient_targets_constructed'])
        self.assertEqual(sum(f['totals']['nonsingular_equations'] for f in report['fields']), 16140)
        first = report['fields'][0]
        self.assertEqual(first['totals']['nonsingular_equations'], 125**2-125)
        self.assertEqual(first['totals']['singular_equations_skipped'], 125)
        for field in report['fields']:
            self.assertEqual(sum(c['equations'] for c in field['module_classes']),
                             field['totals']['nonsingular_equations'])
        candidates = report['fields'][1]['prime_order_extension_j_candidates']
        self.assertEqual(len(candidates), 2)
        for row in candidates:
            self.assertFalse(row['j_in_prime_field'])
            self.assertEqual(row['rational_2_kernels'], 0)
            self.assertEqual(row['mixed_3_kernel_candidates'], 12)
            self.assertEqual(factors(row['group_order']), [[row['group_order'], 1]])

    def test_census_small_sample_is_deterministic(self):
        first, second = census(5, 12), census(5, 12)
        self.assertEqual(first, second)
        self.assertEqual(first['totals']['nonsingular_equations'], 12)
        self.assertEqual(len(first['ordered_curve_data_sha256']), 64)

    def test_prym_saved_batch_reproduces(self):
        report = run_pryms()
        self.assertEqual(json.loads(json.dumps(report)), saved('prym-cover-experiments'))
        self.assertEqual(report['case_count'], 32)
        self.assertEqual(report['cases_with_odd_Prym_prime_absent_from_elliptic_factor'], 23)
        self.assertFalse(report['scope']['Jacobian_DLP_solved'])

    def test_prym_example_counts_with_independent_field_and_point_enumeration(self):
        row = build_case(7, [1, 0, 0, 6, 1])
        self.assertEqual(row['jacobian_orders'], {'C': 1104, 'E': 12, 'D': 92})
        self.assertEqual(row['odd_primes_in_D_order_absent_from_E_order'], [23])
        k = Fp2(7)
        for rhs, key in ((row['C_rhs_ascending'], 'C_genus3'),
                         (row['f_ascending'], 'E_genus1'),
                         (row['D_rhs_ascending'], 'D_genus2')):
            affine = sum(k.mul(y, y) == evaluate(k, rhs, x)
                         for x, y in itertools.product(range(49), repeat=2))
            n = affine+(1 if (len(rhs)-1) % 2 else 2)
            self.assertEqual(n, row['point_counts_over_p_p2_p3'][key][1])
        with self.assertRaises(ValueError):
            build_case(7, [0, 0, 0, 0, 1])

    def test_torus_group_and_scalar_transfer_and_smooth_negative_control(self):
        report = run_tori()
        self.assertEqual(json.loads(json.dumps(report)), saved('torus-transfer-experiments'))
        for row in report['experiments']:
            self.assertEqual(row['group_law_pairs_checked'], (row['p']+1)**2)
            dlp = row['toy_DLP']
            self.assertEqual(dlp['secret'], dlp['recovered_by_torus_enumeration'])
            counter = row['smooth_elliptic_negative_control']['explicit_counterexample']
            self.assertNotEqual(counter['map_P_plus_Q'], counter['map_P_times_map_Q'])
        self.assertFalse(report['scope']['smooth_elliptic_to_torus_reduction_established'])

    def test_torus_cli_and_prym_cli(self):
        for name, result in (('torus_transfer_experiments', 'torus-transfer-experiments'),
                             ('prym_cover_experiments', 'prym-cover-experiments')):
            completed = subprocess.run([sys.executable, '-m', 'research.'+name],
                                       cwd=ROOT, capture_output=True, text=True, check=True)
            self.assertEqual(json.loads(completed.stdout), saved(result))


if __name__ == '__main__':
    unittest.main()
