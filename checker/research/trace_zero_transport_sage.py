"""Evaluate the toy kernel-810 subgroup map using enumerative divisor inversion.

This attempts a real quotient-section/Jacobian match, not an equal-order group
identification. The inverse table remains O(r), and this is not a global
algebraic quotient-map formula or an efficient inversion algorithm.
"""

import argparse
import json
import sys
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from sage.all import GF, PolynomialRing
from sage.env import SAGE_VERSION

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from research.eight_geometry_sage import quartic_curve
from research.quotient_divisor_verification import normalize, verify_source

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")

def run(limit):
    root = Path(__file__).resolve().parents[1]
    path = root / "research/results/eight-trace-zero-inputs.json"
    data = json.loads(path.read_text())
    row = data["sources"][0]
    start = perf_counter()
    ctx = verify_source(data, row, context_only=True, trace_zero_generator=True)
    setup = perf_counter() - start
    K, G = ctx["field"], ctx["group"]
    s = K.gen()

    def elt(cs):
        return sum(K(c) * s**i for i, c in enumerate(cs))

    source_images = []
    for encoded in row["source_seven_coordinate_images_hex"]:
        cs = bytes.fromhex(encoded)
        source_images.append(
            normalize([elt(cs[18 * i : 18 * i + 18]) for i in range(7)])
        )
    assert len(set(source_images)) == 187
    inverse, points, timings = {}, [], []
    current = G.zero()
    for n in range(limit):
        before = perf_counter()
        image = ctx["jacobian_image"](current)
        assert image in source_images and image not in inverse
        inverse[image] = [n] if n == 0 else [n, 373 - n]
        points.append(current)
        timings.append(perf_counter() - before)
        if n % 5 == 0:
            print(
                f"trace-zero inversion {n}/{limit}: {timings[-1]:.3f}s",
                file=sys.stderr,
                flush=True,
            )
        current += ctx["generator"]
    result = {
        "experiment": 2,
        "sage_version": SAGE_VERSION,
        "input_sha256": sha256(path.read_bytes()).hexdigest(),
        "setup_seconds": setup,
        "classes_modulo_sign_evaluated": limit,
        "evaluation_seconds": sum(timings),
        "per_class_seconds": timings,
        "inverse_table_size": len(inverse),
        "inversion_method": "enumeration of the actual trace-zero subgroup modulo sign",
        "inversion_cost_linear_in_r": True,
        "scope": {
            "global_algebraic_quotient_map_constructed": False,
            "scalable_inversion_implemented": False,
            "trace_zero_relation_generation_implemented": False,
            "practical_speedup_established": False,
        },
    }
    if limit < 187:
        result["status"] = "partial_exact_section_matching_control"
        return result
    anchor = min(inverse[source_images[1]])
    image_by_scalar = [None] * 373
    for image, candidates in inverse.items():
        for scalar in candidates:
            image_by_scalar[scalar] = image
    lifts = []
    for n, image in enumerate(source_images):
        candidates = [
            j
            for j in inverse[image]
            if image_by_scalar[(j + anchor) % 373] == source_images[(n + 1) % 373]
        ]
        assert len(candidates) == 1
        lifts.append(candidates[0])
    assert all(j == n * anchor % 373 for n, j in enumerate(lifts))
    # Independently verify all generator additions in the small F49 group.
    F7 = GF(7)
    L = GF(49, "w", modulus=PolynomialRing(F7, "t")([1, 0, 1]))
    C, _, chart = quartic_curve(810, L)
    a, b, c = chart

    def moved(P):
        xx, yy, zz = P
        return C([xx - a * (yy - c * zz) - b * zz, yy - c * zz, zz])

    P = moved([L(0), 1 + 3 * L.gen(), L(1)])
    Q = moved([L(0), 1 + 3 * L.gen() ** 7, L(1)])
    O = next(R for R in C.rational_points() if all(v**7 == v for v in R)).place()
    H = C.jacobian(model="hess", base_div=3 * O.divisor()).group()
    generator = H.point(P.place() - Q.place())
    assert generator != H.zero() and 373 * generator == H.zero()
    # The K embedding may choose w with the other sign. This changes only the
    # permitted global sign of the homomorphism, not the group identities.
    table, current = [], H.zero()
    for _ in range(373):
        table.append(current)
        current += generator
    assert current == H.zero()
    for n in range(373):
        assert table[lifts[n]] + table[anchor] == table[lifts[(n + 1) % 373]]
    assert 37 * table[anchor] == table[lifts[37]]
    result.update(
        {
            "status": "all_373_source_inputs_lifted_and_verified",
            "source_points_lifted": 373,
            "generator_addition_checks": 373,
            "source_generator_image_scalar_in_trace_zero_generator": anchor,
            "source_generator_image_divisor": str(table[anchor].divisor()),
            "Frobenius_sign": -1,
            "Frobenius_reason": "generator=[P-Frob7(P)], all images integer multiples; Frob7^2(P)=P",
            "scalar37_control": True,
            "scope": {
                **result["scope"],
                "actual_subgroup_divisor_transport_verified": True,
            },
        }
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=187, choices=range(1, 188))
    args = parser.parse_args()
    result = run(args.limit)
    result["source_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
