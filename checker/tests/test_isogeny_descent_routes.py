import json
import unittest
from pathlib import Path

from research.isogeny_descent_routes import evaluate_isogeny, on_curve, quotient, run
from research.small_fields import SmallField, ec_add, ec_mul, points


class IsogenyDescentRoutesTests(unittest.TestCase):
    def test_rational_formulas_against_kernel_point_sums_in_other_characteristics(self):
        # An independent oracle: adjoin the missing kernel y-coordinate and
        # compute Velu's point sums using the existing elliptic group law.
        for p in (5, 7, 11):
            k = SmallField(p)
            squares = {k.pow(x, 2) for x in range(p)}
            checked = 0
            for a in range(p):
                for b in range(p):
                    if (4 * a**3 + 27 * b**2) % p == 0:
                        continue
                    for u in range(p):
                        value = (u**3 + a * u + b) % p
                        if value in squares:
                            continue
                        try:
                            aa, bb = quotient(k, a, b, u)
                        except ValueError:
                            continue
                        L = SmallField(p, 2, [(-value) % p, 0, 1])
                        Q, negative_Q = (u, p), (u, L.neg(p))
                        self.assertIsNone(ec_mul(L, 0, a, Q, 3))
                        for P in points(k, [b, a, 0, 1]):
                            plus, minus = (
                                ec_add(L, 0, a, P, Q),
                                ec_add(L, 0, a, P, negative_Q),
                            )
                            xx = L.sub(
                                L.add(P[0], L.add(plus[0], minus[0])), L.mul(2, u)
                            )
                            yy = L.add(P[1], L.add(plus[1], minus[1]))
                            native = evaluate_isogeny(k, a, b, u, P)
                            self.assertEqual(native, (xx, yy))
                            self.assertTrue(on_curve(k, aa, bb, native))
                        checked += 1
                        if checked == 4:
                            break
                    if checked == 4:
                        break
                if checked == 4:
                    break
            self.assertEqual(checked, 4)

    def test_entire_graph_and_saved_report(self):
        root = Path(__file__).resolve().parents[1]
        result = json.loads(json.dumps(run()))
        saved = json.loads(
            (root / "research/results/isogeny-descent-routes.json").read_text()
        )
        self.assertEqual(result, saved)
        summary = result["summary"]
        self.assertEqual(summary["j_vertices"], 9)
        self.assertEqual(summary["directed_edges"], 18)
        self.assertEqual(summary["vertices_by_mixed_kernel_count"], {"0": 6, "12": 3})
        self.assertEqual(summary["point_images_verified"], 6714)
        self.assertEqual(summary["generator_additions_verified"], 6714)
        self.assertEqual(len(summary["j_frobenius_orbits"]), 3)
        self.assertEqual(len(result["routes_to_prior_descent_example"]), 6)
        pairs = {(e["source_j"], e["target_j"]) for e in result["edges"]}
        self.assertTrue(all((v, u) in pairs for u, v in pairs))
        self.assertFalse(result["scope"]["new_Jacobian_maps_evaluated"])
        self.assertFalse(result["scope"]["practical_speedup_established"])

    def test_kernel_identity_exception_and_invalid_root(self):
        k = SmallField(7, 3, [1, 0, 1, 1])
        self.assertIsNone(evaluate_isogeny(k, 248, 42, 19, None))
        with self.assertRaises(ValueError):
            quotient(k, 248, 42, 0)

    def test_trace_zero_report_links_to_the_obstructed_target(self):
        # Sage arithmetic is a separate reproducible experiment. This check
        # guards the boundary: a trace-zero witness is not a quotient map.
        root = Path(__file__).resolve().parents[1]
        report = json.loads(
            (root / "research/results/descent-leads-sage.json").read_text()
        )
        target = report["opposite_sign_target"]
        self.assertEqual(
            target["F7_Jacobian_order_from_full_Frobenius_polynomial"], 315
        )
        self.assertEqual(target["trace_zero_order_from_full_Frobenius_polynomial"], 373)
        self.assertFalse(target["explicit_source_to_this_target_map_evaluated"])
        self.assertFalse(report["scope"]["trace_zero_relation_algorithm_implemented"])
        self.assertEqual(
            report["independent_isogeny_checks"]["rational_point_images"], 6714
        )


if __name__ == "__main__":
    unittest.main()
