"""Independent Sage checks for isogeny routing and a quadratic trace-zero target.

Run with Sage 10.6, from checker/:
sage -python research/descent_leads_sage.py --output /tmp/descent-leads-sage.json
No quotient-to-Jacobian point map is constructed by this experiment.
"""

import argparse
import itertools
import json
from hashlib import sha256
from pathlib import Path

from sage.all import (
    GF,
    ZZ,
    Curve,
    EllipticCurve,
    PolynomialRing,
    hilbert_class_polynomial,
)
from sage.env import SAGE_VERSION

if not __debug__:
    raise RuntimeError("The research certificates require assertions; omit -O")


def run(root):
    route_path = root / "research/results/isogeny-descent-routes.json"
    quotient_path = root / "research/results/quotient-reconstruction.json"
    routes = json.loads(route_path.read_text())
    previous = json.loads(quotient_path.read_text())
    Fp = GF(7)
    R = PolynomialRing(Fp, "t")
    K = GF(343, "z", modulus=R([1, 0, 1, 1]))
    z = K.gen()

    def elt(n):
        return sum(K((n // 7**i) % 7) * z**i for i in range(3))

    nodes = {row["j"]: row for row in routes["nodes"]}
    edges_checked = points_checked = 0
    for edge in routes["edges"]:
        row, target = nodes[edge["source_j"]], nodes[edge["target_j"]]
        E = EllipticCurve(K, [elt(row["a"]), elt(row["b"])])
        ring = PolynomialRing(K, "x")
        phi = E.isogeny(ring.gen() - elt(edge["kernel_x"]))
        assert phi.degree() == 3
        assert tuple(phi.codomain().a_invariants()[3:]) == (
            elt(edge["normalized_codomain_a"]),
            elt(edge["normalized_codomain_b"]),
        )
        target_curve = EllipticCurve(K, [elt(target["a"]), elt(target["b"])])
        scale = elt(edge["codomain_scaling"])
        source_points = sorted(E.points())
        images = []
        for point in source_points:
            image = phi(point)
            if image.is_zero():
                images.append(None)
            else:
                moved = target_curve(scale**2 * image[0], scale**3 * image[1])
                images.append(
                    tuple(int(v.polynomial()[i]) for v in moved[:2] for i in range(3))
                )
        assert len(images) == len(set(images)) == 373
        # Compare native formulas pointwise, independently of point ordering.
        a, b, u = elt(row["a"]), elt(row["b"]), elt(edge["kernel_x"])
        vv, ww = 6 * u * u + 2 * a, 4 * (u**3 + a * u + b)
        for point in source_points:
            if point.is_zero():
                continue
            x, y = point[:2]
            xx = x + vv / (x - u) + ww / (x - u) ** 2
            yy = y * (1 - vv / (x - u) ** 2 - 2 * ww / (x - u) ** 3)
            assert phi(point) == phi.codomain()(xx, yy)
        assert E.cardinality() == target_curve.cardinality() == 373
        edges_checked += 1
        points_checked += len(source_points)

    # Independent CM certificate: all observed j's are roots of the indicated
    # class polynomials; their degrees also account for all nine vertices.
    cm_checks = []
    for disc in (-59, -531):
        H = hilbert_class_polynomial(disc)
        matches = [
            row
            for row in routes["nodes"]
            if row["endomorphism_discriminant_inferred"] == disc
        ]
        assert len(matches) == H.degree()
        assert all(H.change_ring(K)(elt(row["j"])) == 0 for row in matches)
        cm_checks.append(
            {
                "discriminant": disc,
                "class_number": int(H.degree()),
                "hilbert_class_polynomial_coefficients_ascending": list(map(int, H)),
                "matching_j_vertices": sorted(row["j"] for row in matches),
            }
        )

    target = next(
        row
        for row in previous["sources"][0]["targets"]
        if row["kernel_enumeration_index"] == 810
    )
    plane = PolynomialRing(Fp, names=("X", "Y", "Z"))
    X, Y, Z = plane.gens()
    original = sum(
        Fp(c) * X ** m[0] * Y ** m[1] * Z ** m[2]
        for m, c in target["quartic_coefficients_ascending_monomials"]
    )
    for aa, bb, cc in itertools.product(range(7), repeat=3):
        moved = original(X + aa * Y + bb * Z, Y + cc * Z, Z)
        binary = R([moved.monomial_coefficient(Y**i * Z ** (4 - i)) for i in range(5)])
        if (
            moved(0, 0, 1) != 0
            and binary.degree() == 4
            and binary.gcd(binary.derivative()).degree() == 0
        ):
            break
    else:
        raise ArithmeticError("No supported monic projective chart")
    C0 = Curve(moved)
    assert C0.genus() == 3 and C0.is_smooth()
    L = GF(49, "w", modulus=R([1, 0, 1]))
    C = C0.change_ring(L)
    rational = C.rational_points()
    assert len(C0.rational_points()) == 8 and len(rational) == 50
    O = next(point for point in rational if all(v**7 == v for v in point)).place()
    J = C.jacobian(model="hess", base_div=3 * O.divisor())
    group = J.group()
    nonrational = [P for P in rational if any(v**7 != v for v in P)]
    pairs = []
    for point in nonrational:
        conjugate = C([v**7 for v in point])
        D = group.point(point.place() - conjugate.place())
        assert D != group.zero()
        pairs.append((point, conjugate, D))
    P, conjugate, D = pairs[0]
    assert 373 * D == group.zero()
    sigma_D = group.point(conjugate.place() - P.place())
    assert D + sigma_D == group.zero()
    assert 7 * D != group.zero()  # nontrivial independent scalar control
    # P(T)=T^6-29*T^3+343, so det(1-F^2)=P(1)*P(-1).
    polynomial = PolynomialRing(ZZ, "T")(target["frobenius_polynomial_ascending"])
    invariant, anti_invariant = int(polynomial(1)), int(polynomial(-1))
    assert invariant == 315 and anti_invariant == 373
    assert invariant * anti_invariant == 117495

    def coordinates(point):
        return [[int(v.polynomial()[i]) for i in range(2)] for v in point]

    return {
        "scope": {
            "toy_only": True,
            "Sage_version_used": SAGE_VERSION,
            "explicit_elliptic_isogeny_transport_checked": True,
            "trace_zero_Jacobian_divisor_witness_checked": True,
            "quotient_map_to_opposite_sign_target_constructed": False,
            "trace_zero_relation_algorithm_implemented": False,
            "practical_speedup_established": False,
            "novel_vulnerability_established": False,
        },
        "input_sha256": {
            str(p.relative_to(root)): sha256(p.read_bytes()).hexdigest()
            for p in (route_path, quotient_path)
        },
        "independent_isogeny_checks": {
            "directed_edges": edges_checked,
            "rational_point_images": points_checked,
        },
        "CM_class_polynomial_checks": cm_checks,
        "opposite_sign_target": {
            "source_kernel_index": 810,
            "chart_parameters": [aa, bb, cc],
            "F49_modulus_ascending": [1, 0, 1],
            "F7_curve_points": 8,
            "F49_curve_points": 50,
            "F7_Jacobian_order_from_full_Frobenius_polynomial": invariant,
            "F49_Jacobian_order_from_full_Frobenius_polynomial": invariant
            * anti_invariant,
            "trace_zero_order_from_full_Frobenius_polynomial": anti_invariant,
            "nonrational_curve_points_giving_nonzero_trace_zero_divisors": len(pairs),
            "witness_P": coordinates(P),
            "witness_Frobenius_P": coordinates(conjugate),
            "witness_divisor": "[P-Frobenius_7(P)]",
            "witness_order": 373,
            "witness_order_checked_by_Jacobian_scalar_multiplication": True,
            "Frobenius_witness_equals_negative_witness": True,
            "witness_class_modulo_sign_is_F7_rational": True,
            "rational_order_373_subgroup_over_F7": False,
            "explicit_source_to_this_target_map_evaluated": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    result = run(root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
