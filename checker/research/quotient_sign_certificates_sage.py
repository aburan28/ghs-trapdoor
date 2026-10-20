"""Sage: conditional +/- twist certificates from base-field Jacobian points.

The theta construction supplies a geometric quotient and the two possible
Frobenius polynomials. These certificates distinguish those two possibilities;
they are not a general isogeny-class test for arbitrary equal-order Jacobians.
"""

import argparse
import itertools
import json
from hashlib import sha256
from pathlib import Path

from sage.all import GF, Curve, PolynomialRing
from sage.env import SAGE_VERSION

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def run(path):
    source = json.loads(path.read_text())["sources"][0]
    Fp = GF(7)
    plane = PolynomialRing(Fp, names=("X", "Y", "Z"))
    X, Y, Z = plane.gens()
    ring = PolynomialRing(Fp, "t")
    certificates = []
    for row in source["targets"]:
        f = sum(
            Fp(c) * X ** m[0] * Y ** m[1] * Z ** m[2]
            for m, c in row["quartic_coefficients_ascending_monomials"]
        )
        for a, b, c in itertools.product(range(7), repeat=3):
            moved = f(X + a * Y + b * Z, Y + c * Z, Z)
            binary = ring(
                [moved.monomial_coefficient(Y**i * Z ** (4 - i)) for i in range(5)]
            )
            if (
                moved(0, 0, 1) != 0
                and binary.degree() == 4
                and binary.gcd(binary.derivative()).degree() == 0
            ):
                break
        else:
            raise ArithmeticError("No supported monic chart")
        C = Curve(moved)
        points = C.rational_points()
        O, P = points[:2]
        group = C.jacobian(model="hess", base_div=3 * O.place().divisor()).group()
        D = group.point(P.place() - O.place())
        assert D != group.zero()
        rD = 373 * D
        matches = rD == group.zero()
        assert matches == row["matches_Weil_restriction_Frobenius"]
        assert (315 * D == group.zero()) == (not matches)
        certificates.append(
            {
                "kernel_enumeration_index": row["kernel_enumeration_index"],
                "chart_parameters": [a, b, c],
                "point_P": list(map(int, P)),
                "point_O": list(map(int, O)),
                "nonzero_divisor_class": "[P-O]",
                "multiplication_by_373_is_zero": matches,
                "multiplication_by_315_is_zero": not matches,
                "quotient_to_model_Frobenius_lift_sign_under_twist_dichotomy": 1
                if matches
                else -1,
                "sign_matches_prior_full_Frobenius_polynomial": True,
                "extension_field_curve_counts_used_to_form_certificate": False,
                "explicit_quotient_map_Frobenius_conjugation_evaluated": False,
            }
        )
    return {
        "scope": {
            "toy_only": True,
            "Sage_version_used": SAGE_VERSION,
            "conditional_on_the_proven_geometric_quotient_and_twist_dichotomy": True,
            "not_a_general_equal_order_isogeny_test": True,
            "source_models_deduplicated_to_first_galois_orbit": True,
            "quotient_to_Jacobian_map_evaluated_here": False,
            "practical_speedup_established": False,
        },
        "input_sha256": sha256(path.read_bytes()).hexdigest(),
        "summary": {
            "kernel_constructions": len(certificates),
            "positive_lift_sign": sum(
                c["multiplication_by_373_is_zero"] for c in certificates
            ),
            "negative_lift_sign": sum(
                not c["multiplication_by_373_is_zero"] for c in certificates
            ),
        },
        "certificates": certificates,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).with_name("results") / "quotient-reconstruction.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
