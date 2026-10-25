import itertools
import json
import subprocess
import sys
import unittest
from pathlib import Path

from ghs_checker.field import BinaryField
from research.dimension_four_kernels import (
    contains,
    elliptic_form,
    genus_two_form,
    rref_four_planes,
    transform,
)
from research.isogeny_descent_routes import evaluate_isogeny, quotient
from research.small_fields import SmallField, ec_mul, integer_poly_mul
from research.wild_artin_schreier import class_span, span

RESULTS = Path(__file__).resolve().parents[1] / "research/results"


def report(name):
    return json.loads((RESULTS / f"eight-{name}.json").read_text())


class EightExperimentsTests(unittest.TestCase):
    def test_every_saved_route_against_independent_native_formulas(self):
        routes_checked = 0
        for field in report("isogeny-search")["fields"]:
            k = SmallField(field["p"], 3, field["field_modulus_ascending"])
            if k.size <= 343:
                k.tables()
            for row in field["routes"]:
                a, b, u = row["source_a"], row["source_b"], row["kernel_x"]
                self.assertEqual(
                    quotient(k, a, b, u), (row["target_a"], row["target_b"])
                )
                generator = tuple(row["generator"])
                self.assertEqual(
                    evaluate_isogeny(k, a, b, u, generator), tuple(row["image"])
                )
                self.assertIsNone(ec_mul(k, 0, a, generator, row["subgroup_order"]))
                self.assertNotEqual(row["subgroup_order"] % 3, 0)
                routes_checked += 1
        self.assertEqual(routes_checked, 558)

    def test_lagrangian_enumeration_and_recorded_kernels(self):
        counts = [0, 0]
        total = 0
        for rows in rref_four_planes():
            total += 1
            for i, form in enumerate((elliptic_form, genus_two_form)):
                counts[i] += all(
                    (form(a) & b).bit_count() % 2 == 0
                    for a, b in itertools.combinations(rows, 2)
                )
        self.assertEqual((total, counts), (200787, [2295, 2295]))
        for module in report("dimension-four-kernels")["modules"]:
            for kernel in module["kernels"]:
                rows = kernel["basis_rows_binary"]
                for bits in itertools.product((0, 1), repeat=4):
                    v = 0
                    for bit, row in zip(bits, rows):
                        if bit:
                            v ^= row
                    self.assertTrue(
                        contains(transform(v, module["Frobenius_rows_binary"]), rows)
                    )

    def test_artin_schreier_constant_and_pole_normalization(self):
        for n, modulus in ((4, 19), (6, 67)):
            k = BinaryField(n, modulus)
            wp = {k.square(x) ^ x for x in range(k.size)}
            self.assertEqual(wp, {x for x in range(k.size) if k.absolute_trace(x) == 0})
            for c in range(k.size):
                root = k.pow(c, k.size // 2)
                self.assertEqual(k.square(root), c)
                for x in range(k.size):
                    h = k.mul(root, x)
                    self.assertEqual(
                        k.mul(c, k.square(x)) ^ k.mul(root, x), k.square(h) ^ h
                    )

    def test_frobenius_spans_and_conductor_genera(self):
        saved = report("wild-classes")
        self.assertEqual(
            sum(f["families"][0]["reduced_classes_examined"] for f in saved["fields"]),
            8704,
        )
        for field in saved["fields"]:
            n = field["degree"]
            k = BinaryField(n, field["modulus_binary"])
            mask = k.size - 1
            for family in field["families"]:
                self.assertEqual(
                    sum(r["classes"] for r in family["distribution"]),
                    2 * k.size * k.size,
                )
                for row in family["representatives"]:
                    basis = class_span(k, row["a"], row["b"], row["constant_trace_bit"])
                    values = set(span(basis))
                    self.assertEqual(len(values), 2 ** row["span_rank"])
                    for v in values:
                        image = (
                            k.square(v & mask)
                            | (k.square((v >> n) & mask) << n)
                            | (v & (1 << (2 * n)))
                        )
                        self.assertIn(image, values)
                    if row["geometrically_connected"]:
                        expected = sum(
                            bool(v & mask)
                            and (family["family"] == "cubic" or bool((v >> n) & mask))
                            for v in values
                        )
                        self.assertEqual(
                            row["compositum_genus_from_conductors"], expected
                        )

    def test_product_features_do_not_separate_serre_sign_controls(self):
        data = report("obstruction-filter")
        self.assertEqual(len(data["constructions"]), 24)
        self.assertEqual(
            data["theta_product_feature_groups"],
            {"(9, 1, 9)": {"matching": 22, "opposite": 2}},
        )
        self.assertFalse(data["scope"]["normalized_chi18_implemented"])
        self.assertTrue(data["differential_twist_control"]["factor_is_nonsquare"])

    def test_full_prym_and_wild_frobenius_factors(self):
        for case in report("prym")["cases"]:
            p = case["full_Frobenius_polynomials_ascending"]
            self.assertEqual(integer_poly_mul(p["B"], p["D"]), p["C"])
            self.assertEqual(len(p["C"]), 19)
            for n, (c, b, d) in enumerate(
                zip(
                    case["counts_C_over_p_through_p9"],
                    case["counts_B_over_p_through_p9"],
                    case["counts_D_over_p_through_p9"],
                ),
                1,
            ):
                self.assertEqual(c, b + d - 3**n - 1)
        wild = report("wild-cover")
        product = [1]
        for polynomial in wild["quotient_Frobenius_polynomials_ascending"]:
            product = integer_poly_mul(product, polynomial)
        self.assertEqual(product, wild["full_Frobenius_polynomial_H_ascending"])

    def test_actual_subgroup_certificates_keep_their_scope(self):
        genus4 = report("dimension-four-control")
        self.assertEqual(genus4["order17_subgroup_witness"]["order"], 17)
        self.assertTrue(
            genus4["scope"]["actual_genus4_Jacobian_subgroup_image_verified"]
        )
        self.assertFalse(genus4["scope"]["general_Schottky_filter_implemented"])
        prym = report("prym-subgroup")
        self.assertEqual(prym["subgroup_order"], 41)
        self.assertTrue(prym["actual_source_and_image_orders_verified"])
        wild = report("wild-subgroup")
        self.assertTrue(
            wild["actual_nonzero_source_pullback_and_trace_orders_verified"]
        )
        self.assertEqual(wild["trace_image_field"], "F2")
        self.assertFalse(wild["ordinary_elliptic_source"])

    def test_cross_characteristic_sign_split_and_full_transfer(self):
        data = report("cross-characteristic-targets")
        self.assertEqual(data["source_points_transferred_and_additions_checked"], 2137)
        self.assertEqual(
            (data["matching_targets"], data["opposite_sign_targets"]), (3, 9)
        )
        self.assertEqual(len(data["targets"]), 12)
        for target in data["targets"]:
            self.assertEqual(
                target["Jacobian_order"], 2137 if target["matches_source"] else 2259
            )
            self.assertEqual(len(target["full_Frobenius_polynomial_ascending"]), 7)
        self.assertFalse(data["scope"]["divisor_map_to_these_targets_evaluated"])

    def test_trace_transport_and_all_lift_sign_controls(self):
        data = report("trace-zero-transport")
        self.assertEqual(data["source_points_lifted"], 373)
        self.assertEqual(data["generator_addition_checks"], 373)
        self.assertEqual(data["inverse_table_size"], 187)
        self.assertEqual(data["Frobenius_sign"], -1)
        self.assertTrue(data["inversion_cost_linear_in_r"])
        signs = report("lift-signs")
        self.assertEqual(sum(len(row["checks"]) for row in signs["targets"]), 100)
        opposite = next(row for row in signs["targets"] if row["kernel"] == 810)
        nonrational = [check for check in opposite["checks"] if check["sign"] == -1]
        self.assertEqual(len(nonrational), 42)
        self.assertEqual({check["order"] for check in nonrational}, {373})

    def test_wild_comparison_retains_regular_and_irregular_boundaries(self):
        fields = report("wild-classical-comparison")["fields"]
        self.assertEqual(sum(f["regular_comparisons_equal_genus"] for f in fields), 382)
        self.assertEqual(
            sum(
                f["irregular_cases_not_claimed_as_base_field_descents"] for f in fields
            ),
            56,
        )
        for field in fields:
            for row in field["comparisons"]:
                if row["classical_regularity"]:
                    self.assertEqual(
                        row["geometric_span_genus_over_K"], row["classical_GHS_genus"]
                    )
                else:
                    self.assertIsNone(row["classical_GHS_genus"])

    def test_sparse_inversion_failures_are_explicit(self):
        data = report("divisor-inversion-circuit")
        self.assertEqual(len(data["controls"]), 4)
        self.assertEqual(data["variables"], 96)
        for row in data["controls"]:
            self.assertTrue(row["all_forward_equations_verified"])
            self.assertTrue(row["independent_Hess_class_evaluation_matches"])
            self.assertTrue(row["isolated_process"])
        self.assertFalse(data["scalable_inversion_established"])
        self.assertFalse(data["subgroup_enumeration_used"])

    def test_optimized_python_cannot_issue_research_certificates(self):
        for module in (
            "wild_artin_schreier",
            "dimension_four_kernels",
            "wild_classical_comparison",
        ):
            result = subprocess.run(
                [sys.executable, "-O", "-m", "research." + module, "--help"],
                cwd=RESULTS.parents[1],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("omit -O", result.stderr)


if __name__ == "__main__":
    unittest.main()
