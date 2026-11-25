"""Independent small-field checks for the experimental genus-2 audit."""

import itertools
import json
from pathlib import Path
import subprocess
import sys
import unittest

from research.genus2_audit import (
    Fp2, Fp4, ROSENHAIN, BK, audit_bk, audit_rosenhain, branch_screen,
    point, pairings, pair_poly, pmul, peval, richelot, quadratic_roots,
    mapply, madj, mmul, mnormalize, triple_matrix,
)


def affine(f, p):
    return None if p[1] == 0 else f.div(p[0], p[1])


def all_descent_data(k):
    """Independent exhaustive oracle: conjugate Frobenius by every PGL2(K) map."""
    line = [None] + list(range(k.size))
    data = set()
    for triple in itertools.permutations(line, 3):
        g = triple_matrix(k, [point(x) for x in triple])
        m = mmul(k, madj(k, g), tuple(k.sigma(c) for c in g))
        data.add(mnormalize(k, m))
    return data


def brute_branch_condition(k, branches, descent_data):
    s = set(branches)
    for m in descent_data:
        images = {x: affine(k, mapply(k, m, point(None if x is None else k.sigma(x))))
                  for x in branches}
        for omitted in s:
            t = s - {omitted}
            if {images[x] for x in t} == t:
                return True
    return False


def curve_count(f, rhs):
    # Independent square table, rather than calling the square-root implementation.
    squares = {f.mul(y, y) for y in range(f.size)}
    count = 0
    for x in range(f.size):
        value = peval(f, rhs, x)
        count += 1 if value == 0 else 2 if value in squares else 0
    # Smooth projective points at infinity, for degrees 5 and 6.
    count += 1 if len(rhs) == 6 else 2 if rhs[-1] in squares else 0
    return count


def frobenius_polynomial(k, extension, rhs):
    q = k.size
    n1, n2 = curve_count(k, rhs), curve_count(extension, rhs)
    c1 = n1 - q - 1
    c2 = (c1 * c1 + n2 - q*q - 1) // 2
    return (q*q, q*c1, c2, c1, 1)


class Genus2ResearchTests(unittest.TestCase):
    def test_field_arithmetic_roots_and_frobenius(self):
        k = Fp2(7)
        f = Fp4(k)
        squares_k = {k.mul(y, y) for y in range(k.size)}
        for x in range(k.size):
            root = k.sqrt(x)
            self.assertEqual(root is not None, x in squares_k)
            if root is not None:
                self.assertEqual(k.mul(root, root), x)
            lifted = f.sqrt_base(x)
            self.assertEqual(f.mul(lifted, lifted), x)
        for x in range(f.size):
            self.assertEqual(f.sigma(x), f.pow(x, 7))
            if x:
                self.assertEqual(f.mul(x, f.inv(x)), 1)
        self.assertEqual(f.pow(f.q, f.size), f.q)

    def test_branch_screen_against_exhaustive_pgl2_oracle(self):
        k, f = Fp2(3), Fp4(Fp2(3))
        data = all_descent_data(k)
        self.assertEqual(len(data), 3 * (3**2 + 1))
        # All 210 sets of six K-rational branches, with an independent oracle.
        for branches in itertools.combinations([None] + list(range(k.size)), 6):
            actual = bool(branch_screen(f, branches)['branch_condition_candidates'])
            self.assertEqual(actual, brute_branch_condition(k, branches, data))

    def test_known_weak_form_survives_coordinate_change(self):
        k = Fp2(7)
        f = Fp4(k)
        branches = [None, 0, 1, 2, 3, k.elt(1, 1)]
        g = (k.elt(2, 1), k.elt(3, 2), k.elt(1, 3), k.elt(4, 1))
        transformed = [affine(k, mapply(k, g, point(x))) for x in branches]
        for s in (branches, transformed):
            self.assertGreater(branch_screen(f, s)['branch_condition_candidates'], 0)

    def test_nonsplit_branches_against_exhaustive_oracle(self):
        k = Fp2(3)
        f = Fp4(k)
        data = all_descent_data(k)
        # Two rational branches and two conjugate pairs over K; all lie in Fp^4.
        roots = [None, 0]
        polynomials = ([k.neg(f.delta), 0, 1],
                       [k.sub(1, f.delta), k.neg(2), 1])
        for poly in polynomials:
            roots.extend(quadratic_roots(f, poly))
            for x in roots[-2:]:
                self.assertEqual(peval(f, poly, x), 0)
        if len(set(roots)) != 6:
            raise AssertionError("test branches must be distinct")
        oracle = False
        for m in data:
            for omitted in [x for x in roots if x is None or x < k.size]:
                t = set(roots) - {omitted}
                images = {affine(f, mapply(f, m, point(None if x is None else f.sigma(x))))
                          for x in t}
                oracle |= images == t
        self.assertEqual(bool(branch_screen(f, roots)['branch_condition_candidates']), oracle)

    def test_richelot_full_frobenius_polynomial_by_point_enumeration(self):
        k = Fp2(7)
        f = Fp4(k)
        branches = [None, 0, 1, 2, k.elt(1, 1), k.elt(3, 2)]
        rhs = [1]
        for x in branches[1:]:
            rhs = pmul(k, rhs, pair_poly(k, None, x))
        expected = frobenius_polynomial(k, f, rhs)
        checked = 0
        for partition in pairings(tuple(range(6))):
            delta, factors, target = richelot(k, branches, partition)
            if not delta:
                continue
            roots = [x for h in factors for x in quadratic_roots(f, h)]
            self.assertEqual(len(set(roots)), 6)
            for x in roots:
                if x is not None:
                    self.assertEqual(peval(f, target, x), 0)
            self.assertEqual(frobenius_polynomial(k, f, target), expected)
            checked += 1
        self.assertEqual(checked, 15)

    def test_published_bk_obstruction(self):
        for _, p, a in BK:
            report = audit_bk(p, a)
            self.assertTrue(report['quintic_irreducible'])
            self.assertEqual(report['rational_richelot_kernels'], 0)
        # Exhaust the 15 perfect matchings under a 5-cycle, independently.
        cycle = (0, 2, 3, 4, 5, 1)
        for partition in pairings(tuple(range(6))):
            original = {frozenset(pair) for pair in partition}
            image = {frozenset(cycle[x] for x in pair) for pair in partition}
            self.assertNotEqual(original, image)

    def test_saved_report_matches_parameter_audit(self):
        path = Path(__file__).resolve().parents[1] / 'research/results/genus2-audit.json'
        saved = json.loads(path.read_text())
        for (name, p, coefficients), report in zip(ROSENHAIN, saved['generic_and_kummer']):
            current = {'name': name, **audit_rosenhain(p, coefficients)}
            # JSON roundtrip converts tuple partitions to lists.
            self.assertEqual(json.loads(json.dumps(current)), report)
            self.assertEqual(report['nonsingular_genus2_quotients'], 15)
            self.assertEqual(report['source']['branch_condition_candidates'], 0)
            self.assertTrue(all(n['screen']['branch_condition_candidates'] == 0
                                for n in report['neighbors']))

    def test_explicit_genus4_cover_and_frobenius_polynomial(self):
        from research.genus2_cover_example import build_example, F343
        cubic = F343()
        self.assertTrue(all(Fp2.pow(cubic, x, 343) == x for x in range(343)))
        report = build_example()
        self.assertEqual(report['point_counts']['cover_F7_F49_F343_F2401'],
                         [8, 42, 344, 2266])
        self.assertEqual(report['jacobian_group_order'], 2176)
        self.assertTrue(report['cover']['rational_function_identity_verified'])
        self.assertEqual(report['richelot_neighbor']['branch_screen']['branch_condition_candidates'], 0)
        saved_path = Path(__file__).resolve().parents[1] / 'research/results/toy-cover.json'
        self.assertEqual(json.loads(json.dumps(report)), json.loads(saved_path.read_text()))

    def test_custom_input_cli_and_scope(self):
        root = Path(__file__).resolve().parents[1]
        result = subprocess.run([sys.executable, 'research/genus2_audit.py',
                                 '--input', 'research/examples/weak-branch-toy.json'],
                                cwd=root, text=True, capture_output=True, check=True)
        report = json.loads(result.stdout)
        self.assertTrue(report['scope']['branch_test_only'])
        self.assertFalse(report['scope']['descent_map_constructed'])
        self.assertGreater(report['audit']['source']['branch_condition_candidates'], 0)
        self.assertEqual(report['audit']['partitions_tested'], 15)
        with self.assertRaises(ValueError):
            audit_rosenhain(7, [[2, 0], [2, 0], [3, 1]])
        with self.assertRaises(ValueError):
            audit_rosenhain(7, [[2, 0], [3, 0], [8, 1]])
        with self.assertRaises(ValueError):
            audit_bk(7, [0, 0])


if __name__ == '__main__':
    unittest.main()
