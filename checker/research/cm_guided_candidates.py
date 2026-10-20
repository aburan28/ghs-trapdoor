"""Independently recount the bounded Sage CM candidates with native arithmetic."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from .cubic_curve_census import counts_for_module, curve_data
from .genus2_audit import prime64
from .isogeny_descent_routes import evaluate_isogeny, on_curve, quotient
from .small_fields import SmallField, ec_mul


def verify_candidate(row):
    p = row["p"]
    k = SmallField(p, 3, row["modulus_ascending"])
    actual = curve_data(k, row["a"], row["b"])
    if actual is None:
        raise ArithmeticError("Candidate is singular")
    for key in ("a", "b", "j", "group_order", "trace"):
        if actual[key] != row[key]:
            raise ArithmeticError(f"Independent {key} check failed")
    if actual["three_division_roots"] != sorted(row["three_division_roots"]):
        raise ArithmeticError("Independent division polynomial check failed")
    if actual["frobenius_mod_3_representative"] != ((2, 0), (0, 2)):
        raise ArithmeticError("Actual module is not scalar -I")
    if actual["j_in_prime_field"] or not prime64(row["group_order"]):
        raise ArithmeticError("Expected extension-field j and prime order")
    P = tuple(row["subgroup_generator"])
    if k.pow(P[1], 2) != k.add(k.add(k.pow(P[0], 3), k.mul(row["a"], P[0])), row["b"]):
        raise ArithmeticError("Generator is off curve")
    if ec_mul(k, 0, row["a"], P, row["group_order"]) is not None:
        raise ArithmeticError("Generator order check failed")
    orbit = sorted(k.pow(row["j"], p**i) for i in range(3))
    if orbit != row["j_frobenius_orbit_sorted"] or len(set(orbit)) != 3:
        raise ArithmeticError("Frobenius orbit deduplication failed")
    stable2, mixed2 = counts_for_module(p, 2, actual["frobenius_mod_2_representative"])
    stable3, mixed3 = counts_for_module(p, 3, actual["frobenius_mod_3_representative"])
    return {
        **row,
        "independent_field_and_point_count_check": True,
        "independent_division_polynomial_and_generator_checks": True,
        "rational_2_kernels": stable2,
        "mixed_2_kernel_candidates": mixed2,
        "rational_3_kernels": stable3,
        "mixed_3_kernel_candidates": mixed3,
        "quotient_Jacobian_constructed": False,
    }


def run(path=None):
    path = (
        path or Path(__file__).with_name("results") / "cm-guided-candidates-sage.json"
    )
    data = json.loads(path.read_text())
    candidates = [verify_candidate(row) for row in data["candidates"]]
    orbits = {(r["p"], tuple(r["j_frobenius_orbit_sorted"])) for r in candidates}
    classes = {(r["p"], r["trace"]) for r in candidates}
    if len(orbits) != len(candidates):
        raise ArithmeticError("Duplicate Frobenius orbit")
    first = candidates[0]
    k = SmallField(first["p"], 3, first["modulus_ascending"])
    P = tuple(first["subgroup_generator"])
    neighbors = []
    for u in first["three_division_roots"]:
        a, b = quotient(k, first["a"], first["b"], u)
        target = curve_data(k, a, b)
        image = evaluate_isogeny(k, first["a"], first["b"], u, P)
        image_37 = evaluate_isogeny(
            k, first["a"], first["b"], u, ec_mul(k, 0, first["a"], P, 37)
        )
        if (
            not on_curve(k, a, b, image)
            or image is None
            or image_37 != ec_mul(k, 0, a, image, 37)
            or ec_mul(k, 0, a, image, first["group_order"]) is not None
        ):
            raise ArithmeticError("Bounded isogeny transport failed")
        stable, mixed = counts_for_module(
            first["p"], 3, target["frobenius_mod_3_representative"]
        )
        if target["group_order"] != first["group_order"] or mixed != 0:
            raise ArithmeticError(
                "Expected four floor neighbors in this independent class"
            )
        neighbors.append(
            {
                "kernel_x": u,
                "target_a": a,
                "target_b": b,
                "target_j": target["j"],
                "group_order": target["group_order"],
                "rational_3_kernels": stable,
                "mixed_3_kernel_candidates": mixed,
                "source_generator_image": image,
                "scalar_37_relation_verified": True,
                "source_and_image_order": first["group_order"],
                "all_group_points_enumerated_for_this_map": False,
            }
        )
    return {
        "scope": {
            "toy_only": True,
            "point_counts_exhaustive_for_each_selected_equation": True,
            "enumeration_of_all_curves_over_these_fields": False,
            "CM_selection_and_class_polynomials_checked_separately_in_Sage": True,
            "Jacobian_quotient_construction_implemented": False,
            "DLP_solved": False,
            "novel_vulnerability_established": False,
        },
        "sage_input_sha256": sha256(path.read_bytes()).hexdigest(),
        "summary": {
            "verified_curve_models": len(candidates),
            "distinct_j_Frobenius_orbits": len(orbits),
            "elliptic_isogeny_classes_over_stated_fields": len(classes),
            "endomorphism_field_discriminants": sorted(
                {r["endomorphism_field_discriminant"] for r in candidates}
            ),
        },
        "candidates": candidates,
        "first_class_four_isogeny_neighbors": {
            "source_p": first["p"],
            "source_j": first["j"],
            "source_mixed_3_kernel_candidates": 12,
            "neighbors": neighbors,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
