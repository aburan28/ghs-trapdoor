import itertools
import json
from pathlib import Path
import subprocess
import sys
import unittest

from research.weil_quotient_kernels import (
    Cubic523, frobenius_matrix, lagrangians, pairing, poly_gcd_degree,
    poly_pow_mod, rank, rref_planes, row_times, run, screen, verify_published_torsion,
)


def span(rows, ell):
    return {tuple(sum(c[i]*rows[i][j] for i in range(3)) % ell for j in range(6))
            for c in itertools.product(range(ell), repeat=3)}


def exhaustive_membership_oracle(q, ell, matrix):
    f = frobenius_matrix(q, ell, matrix)
    stable = mixed = 0
    for rows, _ in lagrangians(ell, q % ell):
        vectors = span(rows, ell)
        if not all(row_times(v, f, ell) in vectors for v in vectors):
            continue
        stable += 1
        if any(any(v) and sum(any(v[2*i:2*i+2]) for i in range(3)) == 1
               for v in vectors):
            continue
        mixed += 1
    return stable, mixed


class WeilQuotientKernelTests(unittest.TestCase):
    def test_rref_and_lagrangian_enumeration_counts(self):
        self.assertEqual(sum(1 for _ in rref_planes(2)), 1395)
        self.assertEqual(sum(1 for _ in rref_planes(3)), 33880)
        self.assertEqual(len(lagrangians(2, 1)), 135)
        for q in (1, 2):
            self.assertEqual(len(lagrangians(3, q)), 1120)
        # Every F2 isotropic plane from independent triples, without RREF enumeration.
        vectors = [v for v in itertools.product(range(2), repeat=6) if any(v)]
        oracle = set()
        for a, b, c in itertools.combinations(vectors, 3):
            if any((pairing(a, b, 2, 1), pairing(a, c, 2, 1), pairing(b, c, 2, 1))):
                continue
            points = frozenset(span((a, b, c), 2))
            if len(points) == 8:
                oracle.add(points)
        self.assertEqual(oracle, {frozenset(span(rows, 2)) for rows, _ in lagrangians(2, 1)})

    def test_frobenius_pairing_multiplier_and_cubic_power(self):
        basis = [tuple(int(i == j) for j in range(6)) for i in range(6)]
        for q, ell in ((5, 2), (7, 3), (5, 3), (523, 3)):
            for trace in range(ell):
                matrix = [[0, 1], [-pow(q, 3, ell) % ell, trace]]
                f = frobenius_matrix(q, ell, matrix)
                for u, v in itertools.product(basis, repeat=2):
                    self.assertEqual(pairing(row_times(u, f, ell), row_times(v, f, ell), ell, q),
                                     q*pairing(u, v, ell, q) % ell)
                for j, v in enumerate(basis):
                    for _ in range(3):
                        v = row_times(v, f, ell)
                    expected = [0]*6
                    expected[2*(j//2):2*(j//2)+2] = matrix[j % 2]
                    self.assertEqual(v, tuple(expected))

    def test_fast_screen_against_full_span_membership(self):
        for q, ell, matrix in ((7, 2, [[1, 0], [0, 1]]),
                               (7, 3, [[2, 0], [0, 2]]),
                               (5, 3, [[0, 1], [1, 0]])):
            actual = screen(q, ell, matrix)
            self.assertEqual((actual['frobenius_stable_lagrangians'],
                              actual['stable_kernels_passing_no_coordinate_factor_test']),
                             exhaustive_membership_oracle(q, ell, matrix))

    def test_same_trace_and_determinant_different_answers(self):
        scalar = screen(523, 3, [[2, 0], [0, 2]])
        companion = screen(523, 3, [[0, 1], [2, 1]])
        self.assertEqual(scalar['stable_kernels_passing_no_coordinate_factor_test'], 12)
        self.assertEqual(companion['frobenius_stable_lagrangians'], 1)
        self.assertEqual(companion['stable_kernels_passing_no_coordinate_factor_test'], 0)

    def test_irreducible_frobenius_obstruction(self):
        report = screen(7, 3, [[0, 1], [2, 0]])
        self.assertEqual(report['frobenius_stable_lagrangians'], 0)
        report = screen(523, 2, [[0, 1], [1, 1]])
        self.assertEqual(report['frobenius_stable_lagrangians'], 0)

    def test_actual_published_curve_torsion_verification(self):
        verification = verify_published_torsion()
        self.assertTrue(verification['two_division_cubic_irreducible_verified'])
        self.assertTrue(verification['three_division_polynomial_verified'])
        k = Cubic523()
        for a in (1, 523, 524, k.elt(17, 32, 51), k.size-1):
            self.assertEqual(k.mul(a, k.inv(a)), 1)
            self.assertEqual(k.pow(a, k.size), a)
        self.assertEqual(k.pow(523, 3), 521)
        self.assertEqual(poly_gcd_degree(k, [0, 0, 1], [0, 1]), 1)
        self.assertEqual(poly_gcd_degree(k, [1, 0, 1], [0, 1]), 0)
        self.assertEqual(poly_pow_mod(k, [0, 1], 2, [1, 0, 1]), [k.neg(1)])

    def test_subgroup_scope_and_bad_inputs(self):
        report = screen(523, 3, [[2, 0], [0, 2]], subgroup_order=143069629)
        self.assertTrue(report['subgroup_injective_if_claimed_order_is_valid'])
        self.assertFalse(report['scope']['jacobian_over_base_field_verified'])
        self.assertFalse(report['scope']['point_transport_implemented'])
        self.assertFalse(screen(523, 3, [[2, 0], [0, 2]], subgroup_order=3)
                         ['subgroup_injective_if_claimed_order_is_valid'])
        for args in ((3, 3, [[1, 0], [0, 1]]), (7, 5, [[1, 0], [0, 1]]),
                     (7, 3, [[0, 0], [0, 0]]), (7, 3, [[3, 0], [0, 1]])):
            with self.assertRaises(ValueError):
                screen(*args)

    def test_witnesses_and_result_reproducibility(self):
        root = Path(__file__).resolve().parents[1]
        result = run()
        self.assertEqual(json.loads(json.dumps(result)),
                         json.loads((root/'research/results/weil-quotient-kernels.json').read_text()))
        for witness in result['published_example_kernel_screen']['witnesses']:
            self.assertEqual(rank(witness['basis_rows'], 3), 3)
            self.assertEqual(witness['coordinate_factor_intersection_dimensions'], (0, 0, 0))

    def test_module_cli(self):
        root = Path(__file__).resolve().parents[1]
        completed = subprocess.run([sys.executable, '-m', 'research.weil_quotient_kernels'],
                                   cwd=root, text=True, capture_output=True, check=True)
        report = json.loads(completed.stdout)
        self.assertEqual(report['published_example_kernel_screen']['frobenius_stable_lagrangians'], 16)
        self.assertEqual(report['same_curve_two_torsion_obstruction']['frobenius_stable_lagrangians'], 0)


if __name__ == '__main__':
    unittest.main()
