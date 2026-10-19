import importlib.util
import itertools
import json
import subprocess
import sys
import unittest
from pathlib import Path

HAS_FLINT = importlib.util.find_spec("flint") is not None
if HAS_FLINT:
    from research import quotient_reconstruction as qr
    from research.small_fields import SmallField


@unittest.skipUnless(HAS_FLINT, "install research/requirements-reconstruction.txt")
class QuotientReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = qr.run()
        cls.weil = qr.EllipticWeilSet(
            qr.encode(248), qr.encode(42), qr.encode(19), qr.encode(65)
        )
        cls.weils = [cls.weil, cls.weil.conjugate(1), cls.weil.conjugate(2)]
        cls.rows = qr.lagrangians(3, 1)[153][0]

    def test_field_embedding_against_independent_cubic_arithmetic(self):
        k = SmallField(7, 3, [1, 0, 1, 1])
        for a in range(343):
            b = (37 * a + 19) % 343
            self.assertEqual(qr.encode(k.mul(a, b)), qr.encode(a) * qr.encode(b))
            self.assertEqual(qr.encode(k.add(a, b)), qr.encode(a) + qr.encode(b))
        for bad in (-1, 343, True, "3"):
            with self.assertRaises(ValueError):
                qr.encode(bad)
        # fq_default(0) is truthy; checks must use equality, not Python truthiness.
        self.assertTrue(bool(qr.K(0)))
        self.assertEqual(qr.K(0), 0)

    def test_normal_weil_functions_at_independent_translation_probe(self):
        w = self.weil
        P = qr.add(w.A, w.probe, (qr.K(0), qr.encode(176)))
        for ab, cd in itertools.product(w.three, repeat=2):
            total = tuple((a + b) % 3 for a, b in zip(ab, cd))
            cocycle = (
                w.w3(ab, P) * w.w3(cd, qr.add(w.A, P, w.three[ab])) / w.w3(total, P)
            )
            self.assertEqual(
                cocycle, w.zeta ** (2 * (ab[0] * cd[1] - ab[1] * cd[0]) % 3)
            )
        for i, j in itertools.product(range(4), repeat=2):
            cocycle = w.w2(i, P) * w.w2(j, qr.add(w.A, P, w.two[i])) / w.w2(i ^ j, P)
            self.assertEqual(cocycle, -1 if i & 1 and j & 2 else 1)
        for point in w.three.values():
            if point is not None:
                self.assertEqual(
                    tuple(a.frobenius(3) for a in point), qr.negative(point)
                )

    def test_vanishing_projection_seed_is_recovered(self):
        rows = qr.lagrangians(3, 1)[641][0]
        theta = qr.theta_data(self.weils, rows)
        self.assertTrue(all(x == 0 for x in theta))
        seed, theta, gradients, attempts = qr.choose_seed(self.weils, rows)
        self.assertNotEqual(seed, (0, 0, 0))
        self.assertEqual(sum(x != 0 for x in theta), 36)
        self.assertEqual(len(gradients), 28)
        self.assertEqual(attempts[0]["nonzero_theta_nulls"], 0)

    def test_independent_aronhold_systems_reconstruct_same_rational_curve(self):
        _, _, gradients, _ = qr.choose_seed(self.weils, self.rows)
        first = qr.gradient_quartic(gradients)
        second = qr.gradient_quartic(gradients, (63, 38, 52, 41, 13, 26, 19))
        self.assertEqual(first, second)
        self.assertTrue(all(a.frobenius(1) == a for a in first.values()))

    def test_saved_report_reproduces_all_targets_and_subgroup_images(self):
        path = (
            Path(__file__).resolve().parents[1]
            / "research/results/quotient-reconstruction.json"
        )
        self.assertEqual(
            json.loads(json.dumps(self.report)), json.loads(path.read_text())
        )
        for source in self.report["sources"]:
            self.assertEqual(len(source["targets"]), 12)
            self.assertEqual(source["matching_targets"], 11)
            self.assertEqual(source["opposite_sign_targets"], 1)
            mapped = source["selected_kernel_subgroup_map"]
            self.assertEqual(mapped["source_points_evaluated"], 373)
            self.assertEqual(mapped["distinct_projective_images"], 187)
            self.assertEqual(mapped["kernel_translation_checks"], 81)

    def test_opposite_sign_targets_have_exact_subgroup_obstruction(self):
        for source in self.report["sources"]:
            bad = [
                t
                for t in source["targets"]
                if not t["matches_Weil_restriction_Frobenius"]
            ]
            self.assertEqual(len(bad), 1)
            self.assertEqual(bad[0]["point_counts_over_F7_F49_F343"], [8, 50, 257])
            self.assertEqual(bad[0]["jacobian_order"], 315)
            for t in source["targets"]:
                if t["matches_Weil_restriction_Frobenius"]:
                    self.assertEqual(t["point_counts_over_F7_F49_F343"], [8, 50, 431])
                    self.assertEqual(t["jacobian_order"], 373)

    def test_source_hits_are_one_galois_related_restriction_case(self):
        audit = self.report["source_equivalence_audit"]
        self.assertEqual(audit["independent_Weil_restriction_examples"], 1)
        self.assertEqual(audit["j_frobenius_orbit_encodings"], [62, 16, 325])
        self.assertEqual(audit["point_transport_checks"], 373)
        self.assertEqual(len(audit["corresponding_kernel_indices"]), 12)
        self.assertEqual(len({j for i, j in audit["corresponding_kernel_indices"]}), 12)
        k = SmallField(7, 3, [1, 0, 1, 1])
        # Scaling alone is not the asserted curve isomorphism: conjugation matters.
        self.assertNotEqual(k.mul(k.pow(78, 4), 248), 3)

    def test_smoothness_negative_control_and_projective_change(self):
        self.assertFalse(qr.smooth({(4, 0, 0): qr.K(1), (0, 4, 0): qr.K(1)}))
        target = self.report["sources"][0]["targets"][0]
        curve = {
            tuple(m): qr.K(c)
            for m, c in target["quartic_coefficients_ascending_monomials"]
        }
        moved = qr.substitute(
            curve,
            [
                [qr.K(1), qr.K(0), qr.K(0)],
                [qr.K(0), qr.K(1), qr.K(1)],
                [qr.K(0), qr.K(0), qr.K(1)],
            ],
        )
        self.assertTrue(qr.smooth(moved))
        self.assertEqual(
            qr.count_projective_prime(curve, 2), qr.count_projective_prime(moved, 2)
        )

    def test_translation_outside_kernel_changes_the_image(self):
        P = (qr.K(0), qr.encode(176))
        triple = tuple(tuple(x.frobenius(i) for x in P) for i in range(3))
        seed = (0, 0, 0)
        image = qr.normalize(qr.kummer_sections(self.weils, self.rows, seed, triple))
        shifted = (
            qr.add(self.weil.A, triple[0], self.weil.three[1, 0]),
            triple[1],
            triple[2],
        )
        self.assertNotEqual(
            image,
            qr.normalize(qr.kummer_sections(self.weils, self.rows, seed, shifted)),
        )
        negative = tuple(qr.negative(p) for p in triple)
        self.assertEqual(
            image,
            qr.normalize(qr.kummer_sections(self.weils, self.rows, seed, negative)),
        )

    def test_optimized_python_cannot_silently_skip_certificates(self):
        root = Path(__file__).resolve().parents[1]
        result = subprocess.run(
            [sys.executable, "-O", "-m", "research.quotient_reconstruction"],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("omit -O", result.stderr)


if __name__ == "__main__":
    unittest.main()
