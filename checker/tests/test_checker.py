import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import unittest

from ghs_checker import BinaryCurve, BinaryField, analyze_curve, analyze_descent, field_bounds, parse_element, three_isogenies, scan_isogenies
from ghs_checker.cli import load_curve, main
from ghs_checker.field import is_irreducible
from ghs_checker.isogeny import ThreeIsogeny, linear_roots
from ghs_checker.presets import get_preset, list_presets

ROOT = Path(__file__).resolve().parents[1]


def all_points(curve):
    f = curve.field
    return [None] + [(x, y) for x in range(f.size) for y in range(f.size) if curve.contains((x, y))]


def trial_irreducible(poly):
    """Independent trial-division oracle for small binary polynomials."""
    degree = poly.bit_length() - 1
    for d in range(1, degree // 2 + 1):
        for candidate in range(1 << d, 1 << (d + 1)):
            value = poly
            while value and value.bit_length() >= candidate.bit_length():
                value ^= candidate << (value.bit_length() - candidate.bit_length())
            if value == 0:
                return False
    return degree > 0


class FieldTests(unittest.TestCase):
    def test_input_encodings_and_rejection(self):
        for value in (0x83, '0x83', '0b10000011', '131', 'x^7+x+1', 'z^7+z+1'):
            self.assertEqual(parse_element(value), 0x83)
        self.assertEqual(parse_element('x+x+1'), 1)
        for value in ('x^7;print(1)', '-1', '__import__("os")', True):
            with self.assertRaises(ValueError):
                parse_element(value)
        for degree, poly in ((4, 0x15), (3, 0x13), (0, 1)):
            with self.assertRaises(ValueError):
                BinaryField(degree, poly)

    def test_irreducibility_against_trial_division(self):
        for degree in range(2, 8):
            for poly in range((1 << degree) | 1, 1 << (degree + 1), 2):
                self.assertEqual(is_irreducible(poly), trial_irreducible(poly), hex(poly))

    def test_field_arithmetic(self):
        f = BinaryField(5, 0x25)
        for a in range(32):
            self.assertEqual(f.square(a), f.mul(a, a))
            self.assertEqual(f.pow(a, 32), a)
            if a:
                self.assertEqual(f.mul(a, f.inv(a)), 1)
            self.assertIn(f.absolute_trace(a), (0, 1))
            for b in range(32):
                self.assertEqual(f.mul(a, b), f.mul(b, a))
                self.assertEqual(f.square(a ^ b), f.square(a) ^ f.square(b))


class GHSTests(unittest.TestCase):
    def test_published_teske_pair(self):
        for name, expected in [('teske-secret', (4, 8)), ('teske-public', (7, 64))]:
            curve, order, meta = load_curve(get_preset(name))
            result = analyze_descent(curve, 23, order)
            self.assertEqual((result['m'], result['genus']), expected)
            self.assertEqual(result['subgroup_capacity']['status'], 'not_excluded_by_size')
            self.assertFalse(analyze_curve(curve)['scope']['descent_map_constructed'])

    def test_all_nist_presets_and_koblitz_subgroup_obstruction(self):
        for name, preset in list_presets().items():
            if name.startswith('teske-'):
                continue
            curve, order, metadata = load_curve(preset)
            result = analyze_curve(curve, order)
            self.assertEqual(len(result['descents']), 1)
            d = result['descents'][0]
            if name.startswith('K-'):
                self.assertEqual((d['m'], d['genus']), (1, 1))
                self.assertEqual(d['screening'], 'subgroup_size_obstruction')
            else:
                self.assertEqual(d['m'], preset['degree'])
                self.assertEqual(d['screening'], 'above_genus_threshold')
            point = preset['base_point']
            point = parse_element(point['x']), parse_element(point['y'])
            self.assertTrue(curve.contains(point), name)
            self.assertIsNone(curve.multiply(order, point), name)

    def test_published_nist_field_wide_bounds(self):
        for degree, minimum in [(163, 163), (233, 30), (283, 95), (409, 205), (571, 115)]:
            result = field_bounds(degree)
            self.assertEqual(result['minimum_m_gt_one'], minimum)
            self.assertEqual(result['minimum_genus_for_m_gt_one'], (1 << (minimum - 1)) - 1)

    def test_both_genus_branches_and_even_degree_regularity(self):
        f = BinaryField(7, 0x83)
        seen = set()
        for b in range(1, f.size):
            curve = BinaryCurve(f, 0, b)
            result = analyze_descent(curve, 1)
            if result['m'] == 4:
                seen.add(result['genus'])
                # Independent odd-degree genus criterion via relative trace.
                self.assertEqual(result['genus'], 7 if f.relative_trace(b, 1) == 0 else 8)
        self.assertEqual(seen, {7, 8})
        f = BinaryField(4, 0x13)
        trace_one = next(a for a in range(16) if f.absolute_trace(a) == 1)
        result = analyze_descent(BinaryCurve(f, trace_one, 1), 1)
        self.assertEqual(result['screening'], 'regularity_condition_failed')
        self.assertIsNone(result['genus'])

    def test_frobenius_and_basis_independence(self):
        f = BinaryField(4, 0x13)
        alt = BinaryField(4, 0x19)
        # An explicit field isomorphism: map x to a root of x^4+x+1 in alt.
        root = next(t for t in range(16) if alt.pow(t, 4) ^ t ^ 1 == 0)
        def image(value):
            out, power = 0, 1
            for i in range(4):
                if value & (1 << i):
                    out ^= power
                power = alt.mul(power, root)
            return out
        for b in range(1, 16):
            for s in (1, 2):
                d = analyze_descent(BinaryCurve(f, 0, b), s)
                frob = analyze_descent(BinaryCurve(f, 0, f.square(b)), s)
                transported = analyze_descent(BinaryCurve(alt, 0, image(b)), s)
                self.assertEqual((d['m'], d['genus']), (frob['m'], frob['genus']))
                self.assertEqual((d['m'], d['genus']), (transported['m'], transported['genus']))

    def test_subfields_and_scope(self):
        curve, _, _ = load_curve(get_preset('teske-secret'))
        self.assertEqual([d['subfield_degree'] for d in analyze_curve(curve)['descents']], [1, 7, 23])
        self.assertFalse(analyze_curve(curve, subfields=[23])['all_proper_subfields_checked'])
        self.assertEqual(analyze_curve(BinaryCurve(BinaryField(1, 3), 0, 1))['descents'], [])
        with self.assertRaises(ValueError):
            analyze_descent(curve, 3)
        with self.assertRaises(ValueError):
            analyze_descent(curve, 161)


class IsogenyTests(unittest.TestCase):
    def test_root_factorization_against_exhaustive_search(self):
        for degree, modulus in [(2, 7), (3, 11), (4, 19), (5, 37), (6, 67)]:
            f = BinaryField(degree, modulus)
            for b in range(1, f.size):
                poly = [b, 0, 0, 1, 1]
                self.assertEqual(linear_roots(poly, f, False), linear_roots(poly, f, True))

    def test_isogeny_maps_point_counts_and_addition(self):
        count, nonrational_kernels = 0, 0
        for degree, modulus in [(2, 7), (3, 11), (4, 19)]:
            f = BinaryField(degree, modulus)
            for a in range(f.size):
                for b in range(1, f.size):
                    curve = BinaryCurve(f, a, b)
                    for phi in three_isogenies(curve):
                        points = all_points(curve)
                        self.assertEqual(len(points), phi.target.count_points())
                        kernel = [p for p in points if p is not None and p[0] == phi.kernel_x]
                        nonrational_kernels += not bool(kernel)
                        for p in points:
                            self.assertTrue(phi.target.contains(phi.apply(p)))
                        # Every pair over the smallest fields, sampled pairs in GF(16).
                        pairs = ((p, q) for p in points for q in points) if degree <= 3 else zip(points, reversed(points))
                        for p, q in pairs:
                            self.assertEqual(phi.apply(curve.add(p, q)), phi.target.add(phi.apply(p), phi.apply(q)))
                        count += 1
        self.assertGreater(count, 100)
        self.assertGreater(nonrational_kernels, 0)

    def test_search_finds_explicit_lower_genus_path(self):
        data = json.loads((ROOT / 'examples/toy-isogeny.json').read_text())
        curve, order, meta = load_curve(data)
        self.assertEqual(curve.count_points(), 114)
        point = parse_element(data['base_point']['x']), parse_element(data['base_point']['y'])
        self.assertIsNone(curve.multiply(order, point))
        result = scan_isogenies(curve, order, max_curves=8, max_depth=2, metadata=meta)
        candidates = [node for node in result['nodes'] if node['index'] in result['candidate_nodes']]
        self.assertTrue(candidates)
        for node in candidates:
            target, image = curve, point
            for u in node['path_kernel_x']:
                phi = ThreeIsogeny(target, parse_element(u))
                image, target = phi.apply(image), phi.target
            self.assertIsNotNone(image)
            self.assertIsNone(target.multiply(order, image))
            self.assertEqual(target.count_points(), curve.count_points())
            for scalar in range(order):
                source_point, destination_point = curve.multiply(scalar, point), point
                target_step = curve
                for u in node['path_kernel_x']:
                    phi = ThreeIsogeny(target_step, parse_element(u))
                    source_point = phi.apply(source_point)
                    target_step = phi.target
                self.assertEqual(source_point, target.multiply(scalar, image))

    def test_search_bounds_and_incompatible_subgroup(self):
        curve = BinaryCurve(BinaryField(7, 0x83), 1, 0x18)
        result = scan_isogenies(curve, 19, max_curves=1, max_depth=2)
        self.assertEqual(len(result['nodes']), 1)
        self.assertTrue(result['search']['curve_budget_reached'])
        self.assertTrue(result['search']['bounded_search_only'])
        self.assertFalse(result['scope']['security_certification'])
        with self.assertRaises(ValueError):
            scan_isogenies(curve, 3)


class CLITests(unittest.TestCase):
    def invoke(self, args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = main(args)
        return code, out.getvalue(), err.getvalue()

    def test_json_and_console(self):
        code, out, err = self.invoke(['check', '--preset', 'teske-secret', '--subfield-degree', '23', '--json'])
        self.assertEqual(code, 0, err)
        result = json.loads(out)
        self.assertEqual(result['descents'][0]['genus'], 8)
        self.assertFalse(result['scope']['security_certification'])
        code, out, err = self.invoke(['check', '--preset', 'sect163r2'])
        self.assertEqual(code, 0, err)
        self.assertIn('Classical GHS screening', out)
        code, out, err = self.invoke(['scan', '--input', str(ROOT / 'examples/toy-isogeny.json'), '--json'])
        self.assertEqual(code, 0, err)
        self.assertTrue(json.loads(out)['candidate_nodes'])

    def test_error_paths_and_unsupported_inputs(self):
        for args in [
            ['check', '--degree', '4', '--modulus', 'x^4+x^2+1', '--a', '0', '--b', '1'],
            ['check', '--degree', '3', '--modulus', 'x^3+x+1', '--a', '0', '--b', '0'],
            ['check', '--preset', 'B-163', '--b', '1'],
            ['check', '--preset', 'B-163', '--subfield-degree', '3'],
            ['check', '--preset', 'unknown'],
            ['check'],
            ['scan', '--preset', 'B-163', '--max-curves', '0'],
        ]:
            code, out, err = self.invoke(args + ['--json'])
            self.assertEqual(code, 2, args)
            self.assertIn('error', json.loads(err))
        with self.assertRaises(ValueError):
            load_curve({'characteristic': 3})
        with self.assertRaises(ValueError):
            load_curve({'model': 'supersingular'})

    def test_module_entry_point(self):
        result = subprocess.run([sys.executable, '-m', 'ghs_checker', 'field', '--degree', '233', '--json'],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['minimum_m_gt_one'], 30)


if __name__ == '__main__':
    unittest.main()
