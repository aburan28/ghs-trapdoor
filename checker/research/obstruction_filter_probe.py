"""Probe theta-product filters without misidentifying them as normalized chi18.

All 24 saved kernel constructions are recomputed. Predictions use only theta
data; full Frobenius polynomials are consulted afterwards as independent labels.
"""

import argparse
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from .quotient_reconstruction import (
    EllipticWeilSet,
    K,
    choose_seed,
    coefficients,
    encode,
)

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")

def minimal_field(value):
    return next(d for d in (1, 2, 3, 6, 9, 18) if value.frobenius(d) == value)


def run():
    path = Path(__file__).with_name("results") / "quotient-reconstruction.json"
    previous = json.loads(path.read_text())
    rows = []
    start = perf_counter()
    for source in previous["sources"]:
        x1, x2 = source["source_torsion_x_coordinates"]
        w = EllipticWeilSet(
            encode(source["source_a"]),
            encode(source["source_b"]),
            encode(x1),
            encode(x2),
        )
        weils = [w, w.conjugate(1), w.conjugate(2)]
        for target in source["targets"]:
            before = perf_counter()
            seed, theta, gradients, _ = choose_seed(weils, target["kernel_basis_rows"])
            product = K(1)
            for value in theta:
                if value != 0:
                    product *= value / theta[0]
            degree = minimal_field(product)
            character = product ** ((7**degree - 1) // 2)
            assert character in (K(1), K(-1))
            root = product.sqrt()
            row = {
                "source_a": source["source_a"],
                "source_b": source["source_b"],
                "kernel": target["kernel_enumeration_index"],
                "seed": seed,
                "nonzero_theta_values": sum(v != 0 for v in theta),
                "odd_gradients": len(gradients),
                "product_squared_theta_ratios_coefficients": coefficients(product),
                "product_minimal_field_degree": degree,
                "product_square_character_in_minimal_field": 1
                if character == K(1)
                else -1,
                "chosen_root_minimal_field_degree": minimal_field(root),
                "root_Frobenius_multiplier_coefficients": coefficients(
                    root.frobenius(1) / root
                ),
                "evaluation_seconds": perf_counter() - before,
                "normalized_chi18_status": "unresolved_theta_root_and_differential_normalization",
            }
            # These are comparison labels, never inputs to the product calculation.
            row["independent_full_polynomial_matches_source"] = target[
                "matches_Weil_restriction_Frobenius"
            ]
            rows.append(row)
    groups = {}
    for row in rows:
        key = (
            row["product_minimal_field_degree"],
            row["product_square_character_in_minimal_field"],
            row["chosen_root_minimal_field_degree"],
        )
        groups.setdefault(str(key), Counter())[
            row["independent_full_polynomial_matches_source"]
        ] += 1
    delta = 3
    assert pow(delta, 3, 7) == 6
    # A quadratic twist changes a differential volume by (sqrt(delta))^3.
    # Weight 18 therefore contributes delta^(-27), a nonsquare here.
    factor = pow(delta, -27, 7)
    assert pow(factor, 3, 7) == 6
    return {
        "experiment": 5,
        "constructions": rows,
        "theta_product_feature_groups": {
            key: {"matching": c[True], "opposite": c[False]}
            for key, c in groups.items()
        },
        "differential_twist_control": {
            "delta": delta,
            "weight": 18,
            "volume_twist_exponent": 3,
            "missing_factor_delta_power": -27,
            "factor_F7": factor,
            "factor_is_nonsquare": True,
        },
        "input_sha256": sha256(path.read_bytes()).hexdigest(),
        "elapsed_seconds": perf_counter() - start,
        "scope": {
            "all_24_theta_constructions_recomputed": True,
            "normalized_chi18_implemented": False,
            "reliable_early_rejection_filter_certified": False,
            "failed_normalization_is_not_a_curve_exclusion": True,
            "practical_speedup_established": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    result["source_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
