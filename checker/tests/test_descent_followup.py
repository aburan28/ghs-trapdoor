import copy
import itertools
import json
import subprocess
import sys
import unittest
from hashlib import sha256
from math import gcd
from pathlib import Path

from ghs_checker.field import BinaryField
from research.artin_schreier_classes import (
    add_polynomials,
    coboundary,
    evaluate,
    reduce_polynomial,
)
from research.artin_schreier_classes import (
    run as wild_run,
)
from research.cm_guided_candidates import run as cm_run
from research.cm_guided_candidates import verify_candidate
from research.cubic_curve_census import curve_data
from research.small_fields import SmallField

try:
    from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
except ImportError:
    fq_default_ctx = None


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "research/results"


class DescentFollowupTests(unittest.TestCase):
    def test_artin_schreier_classes_against_exhaustive_coboundary_oracle(self):
        field = BinaryField(3, 0b1011)
        boundaries = {
            tuple(coboundary(field, h)) for h in itertools.product(range(8), repeat=2)
        }
        self.assertEqual(len(boundaries), 32)
        fibers = {}
        for f in itertools.product(range(8), repeat=3):
            normal, correction = reduce_polynomial(field, f)
            fibers[tuple(normal)] = fibers.get(tuple(normal), 0) + 1
            self.assertIn(tuple(add_polynomials(f, normal)), boundaries)
            self.assertEqual(reduce_polynomial(field, normal)[0], normal)
            for x in range(8):
                y_shift = evaluate(field, correction, x)
                for y in range(8):
                    source = field.square(y) ^ y == evaluate(field, f, x)
                    Y = y ^ y_shift
                    target = field.square(Y) ^ Y == evaluate(field, normal, x)
                    self.assertEqual(source, target)
        self.assertEqual(len(fibers), 16)
        self.assertEqual(set(fibers.values()), {32})

    def test_wild_cover_report_and_subgroup_obstruction(self):
        actual = wild_run()
        self.assertEqual(
            actual, json.loads((RESULTS / "artin-schreier-classes.json").read_text())
        )
        self.assertEqual(actual["curve_counts_F8_F64_F512"], [9, 65, 513])
        self.assertEqual(actual["F8_Jacobian_order"], 65)
        self.assertEqual(actual["F512_Jacobian_order"], 262145)
        self.assertEqual(actual["trace_zero_factor_order"], 37 * 109)
        self.assertEqual(
            actual["extension_group_odd_primes_absent_from_base_Jacobian_order"],
            [37, 109],
        )
        self.assertTrue(actual["supersingular_genus_2_toy_control"])
        self.assertFalse(actual["scope"]["norm_conorm_subgroup_map_evaluated"])
        self.assertEqual(
            sum(r["affine_point_transports_checked"] for r in actual["experiments"]),
            16 * 512,
        )

    def test_artin_schreier_frobenius_equivariance_and_limits(self):
        data = json.loads((RESULTS / "artin-schreier-classes.json").read_text())
        field = BinaryField(9, int(data["field_F512"]["modulus_hex"], 16))
        for row in data["experiments"]:
            for n in (1, 3, 6):
                conjugate = [
                    field.frobenius(c, n) for c in row["raw_polynomial_ascending"]
                ]
                normal, _ = reduce_polynomial(field, conjugate)
                self.assertEqual(
                    normal,
                    [
                        field.frobenius(c, n)
                        for c in row["reduced_polynomial_ascending"]
                    ],
                )
        with self.assertRaises(ValueError):
            reduce_polynomial(BinaryField(2, 0b111), [1, 1])
        with self.assertRaises(ValueError):
            reduce_polynomial(field, [512])
        optimized = subprocess.run(
            [sys.executable, "-O", "-m", "research.artin_schreier_classes"],
            check=False,
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(optimized.returncode, 0)
        self.assertIn("Research certificates require assertions", optimized.stderr)

    def test_new_curve_candidates_recount_and_deduplicate(self):
        result = json.loads(json.dumps(cm_run()))
        self.assertEqual(
            result, json.loads((RESULTS / "cm-guided-candidates.json").read_text())
        )
        self.assertEqual(result["summary"]["distinct_j_Frobenius_orbits"], 6)
        self.assertEqual(
            result["summary"]["elliptic_isogeny_classes_over_stated_fields"], 3
        )
        self.assertEqual(
            result["summary"]["endomorphism_field_discriminants"], [-563, -211, -107]
        )
        for row in result["candidates"]:
            self.assertEqual(row["mixed_3_kernel_candidates"], 12)
            self.assertEqual(row["rational_2_kernels"], 0)
            self.assertFalse(row["quotient_Jacobian_constructed"])

    def test_candidate_checks_reject_false_counts_and_points(self):
        row = json.loads((RESULTS / "cm-guided-candidates-sage.json").read_text())[
            "candidates"
        ][0]
        wrong = copy.deepcopy(row)
        wrong["trace"] += 9
        with self.assertRaises(ArithmeticError):
            verify_candidate(wrong)
        wrong = copy.deepcopy(row)
        wrong["subgroup_generator"][1] = 0
        with self.assertRaises(ArithmeticError):
            verify_candidate(wrong)

    def test_four_neighbors_have_different_integral_torsion_modules(self):
        data = json.loads((RESULTS / "cm-guided-candidates.json").read_text())
        source = data["candidates"][0]
        k = SmallField(13, 3, source["modulus_ascending"])
        rows = data["first_class_four_isogeny_neighbors"]["neighbors"]
        self.assertEqual(len(rows), 4)
        self.assertEqual(len({r["target_j"] for r in rows}), 4)
        for row in rows:
            actual = curve_data(k, row["target_a"], row["target_b"])
            self.assertEqual(actual["trace"], source["trace"])
            self.assertEqual(actual["group_order"], 2281)
            self.assertEqual(len(actual["three_division_roots"]), 1)
            self.assertEqual(row["mixed_3_kernel_candidates"], 0)
            self.assertTrue(row["scalar_37_relation_verified"])

    def test_sign_certificates_have_the_conditional_scope_and_correct_twist(self):
        data = json.loads(
            (RESULTS / "quotient-sign-certificates-sage.json").read_text()
        )
        previous = RESULTS / "quotient-reconstruction.json"
        self.assertEqual(
            data["input_sha256"], sha256(previous.read_bytes()).hexdigest()
        )
        self.assertEqual(
            data["summary"],
            {
                "kernel_constructions": 12,
                "positive_lift_sign": 11,
                "negative_lift_sign": 1,
            },
        )
        self.assertEqual(gcd(373, 315), 1)
        negative = [
            r for r in data["certificates"] if not r["multiplication_by_373_is_zero"]
        ]
        self.assertEqual([r["kernel_enumeration_index"] for r in negative], [810])
        self.assertTrue(data["scope"]["not_a_general_equal_order_isogeny_test"])

    @unittest.skipUnless(fq_default_ctx is not None, "requires optional python-flint")
    def test_nonenumerative_inverse_with_independent_field_and_polynomial_checks(self):
        report = json.loads(
            (RESULTS / "symmetric-kummer-inversion-sage.json").read_text()
        )
        modulus = report["splitting_field_modulus_ascending"]
        field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(7)(modulus), var="v")
        ring = fq_default_poly_ctx(field)
        image_w = min(
            (r for r, _ in ring([1, 0, 1]).roots()), key=lambda r: tuple(r.to_list())
        )
        embed = lambda v: field(v[0]) + field(v[1]) * image_w
        zeta_basis = report["zeta_basis_on_120_symmetric_products"]
        triples = list(itertools.combinations_with_replacement(range(8), 3))

        def zeta(points):
            B = [
                [field(1), x, y, x * x, x * y, y * y, x**3, x * x * y]
                for x, y in points
            ]
            products = []
            for ijk in triples:
                total = field(0)
                for indices in itertools.permutations(ijk):
                    total += B[0][indices[0]] * B[1][indices[1]] * B[2][indices[2]]
                products.append(total)
            return [
                sum((field(c) * v for c, v in zip(row, products)), field(0))
                for row in zeta_basis
            ]

        self.assertFalse(report["scope"]["Jacobian_group_enumeration_used"])
        self.assertFalse(
            report["scope"][
                "calibration_from_elliptic_quotient_theta_coordinates_implemented"
            ]
        )
        self.assertEqual(
            (
                report["diagonal_matrix_rank"],
                report["diagonal_vanishing_space_dimension"],
            ),
            (112, 8),
        )
        quartic = {tuple(m): embed(c) for m, c in report["transformed_quartic"]}
        t = ring.gen()
        for row in report["experiments"]:
            self.assertEqual(row["zero_dimensional_fiber_degree"], 12)
            self.assertEqual(row["solver_steps"]["maximum_univariate_degree"], 3)
            target = [embed(c) for c in row["input_eight_projective_values_F49"]]
            divisors = [
                [tuple(field(v) for v in P) for P in D]
                for D in row["recovered_divisor_supports_F7_12"]
            ]
            conic = [field(v) for v in row["residual_conic_coefficients_F7_12"]]
            self.assertEqual(len(divisors), 2)
            for divisor in divisors:
                image = zeta(divisor)
                pivot = next(i for i, v in enumerate(target) if v != 0)
                self.assertTrue(
                    all(
                        v * target[pivot] == target[i] * image[pivot]
                        for i, v in enumerate(image)
                    )
                )
                for x, y in divisor:
                    self.assertEqual(
                        sum(
                            (c * x**i * y**j for (i, j, _), c in quartic.items()),
                            field(0),
                        ),
                        field(0),
                    )
                    self.assertEqual(
                        conic[0] + conic[1] * x + conic[2] * y + conic[3] * x * x,
                        field(0),
                    )
            yy = -(conic[0] + conic[1] * t + conic[3] * t * t) / conic[2]
            actual = sum(
                (c * t**i * yy**j for (i, j, _), c in quartic.items()), ring.zero()
            )
            expected = ring.one()
            for divisor in divisors:
                for x, _ in divisor:
                    expected *= t - x
            self.assertEqual(actual / actual.leading_coefficient(), expected)
            repeated = [divisors[0][0], divisors[0][0], divisors[0][1]]
            self.assertTrue(all(v == 0 for v in zeta(repeated)))


if __name__ == "__main__":
    unittest.main()
