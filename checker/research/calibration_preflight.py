"""Freeze a bounded source frame and unlabeled inputs for the calibration design.

This runs only source-side rank/frame checks. It does not fit a calibration,
evaluate target divisors, or execute the proposed transfer experiment.
Requires the existing optional python-flint reconstruction backend.
"""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from .quotient_reconstruction import (
    MODULUS,
    EllipticWeilSet,
    K,
    add,
    choose_seed,
    coefficients,
    encode,
    kummer_sections,
    lagrangians,
    normalize,
    solve,
)
from .small_fields import SmallField, ec_add


def pivots(rows):
    rows = [list(row) for row in rows]
    result, rank = [], 0
    for j in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][j] != 0), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][j]
        rows[rank] = [v / scale for v in rows[rank]]
        for i in range(len(rows)):
            if i == rank:
                continue
            scale = rows[i][j]
            rows[i] = [u - scale * v for u, v in zip(rows[i], rows[rank])]
        result.append(j)
        rank += 1
        if rank == len(rows):
            break
    return result


def run():
    w = EllipticWeilSet(encode(248), encode(42), encode(19), encode(65))
    weils = [w, w.conjugate(1), w.conjugate(2)]
    rows = lagrangians(3, 1)[153][0]
    seed, _theta, _gradients, _attempts = choose_seed(weils, rows)
    G, current = (K(0), encode(176)), None
    sample, points = [], []
    digest = sha256()
    reverse = {encode(n): n for n in range(343)}
    for _ in range(32):
        current = add(w.A, current, G)
        points.append([reverse[v] for v in current])
        triple = tuple(tuple(v.frobenius(i) for v in current) for i in range(3))
        image = normalize(kummer_sections(weils, rows, seed, triple))
        sample.append(image)
        digest.update(bytes(c for v in image for c in coefficients(v)))
    columns = pivots(sample)
    if len(columns) != 8:
        raise ArithmeticError(
            "The frozen 32-point pilot does not have the expected rank eight"
        )
    vectors = [tuple(row[i] for i in columns) for row in sample]
    basis = []
    for i, row in enumerate(vectors):
        if len(pivots([vectors[j] for j in basis] + [row])) > len(basis):
            basis.append(i)
        if len(basis) == 8:
            break
    transpose = list(map(list, zip(*(vectors[i] for i in basis))))
    ninth = None
    weights = None
    for i, row in enumerate(vectors):
        if i in basis:
            continue
        candidate = solve(transpose, list(row))
        if all(c != 0 for c in candidate):
            ninth, weights = i, candidate
            break
    if ninth is None:
        raise ArithmeticError(
            "The fixed pilot did not yield a projective nine-point frame"
        )
    frame = basis + [ninth]
    if any(
        sum(weights[j] * vectors[basis[j]][i] for j in range(8)) != vectors[ninth][i]
        for i in range(8)
    ):
        raise ArithmeticError("Frame-coordinate reconstruction failed")
    for omitted in range(9):
        if len(pivots([vectors[frame[j]] for j in range(9) if j != omitted])) != 8:
            raise ArithmeticError("A frame minor is singular")

    # Freeze coordinate inputs without calculating their discrete-log labels
    # or target images. Reserve both signs of every pilot point and both
    # signs of every challenge and its P+G auxiliary point.
    small = SmallField(7, 3, [1, 0, 1, 1])
    used_x = {P[0] for P in points}
    challenges = []
    counter = 0
    while len(challenges) < 32 and counter < 10000:
        raw = sha256(f"ghs-calibration-holdout-v1:{counter}".encode()).digest()
        counter += 1
        x = int.from_bytes(raw[:4], "big") % 343
        rhs = small.add(small.add(small.pow(x, 3), small.mul(248, x)), 42)
        y = small.pow(rhs, 86)  # 343 == 3 mod 4
        if small.pow(y, 2) != rhs or y == 0:
            continue
        if raw[4] & 1:
            y = small.neg(y)
        P = (x, y)
        auxiliary = ec_add(small, 0, 248, P, (0, 176))
        if (
            auxiliary is None
            or x in used_x
            or auxiliary[0] in used_x
            or auxiliary[0] == x
        ):
            continue
        case_id = sha256(json.dumps(P, separators=(",", ":")).encode()).hexdigest()
        challenges.append({"case_id": case_id, "point": P, "point_plus_G": auxiliary})
        used_x.update((x, auxiliary[0]))
    if len(challenges) != 32:
        raise ArithmeticError("The frozen bounded challenge stream exhausted its limit")
    results = Path(__file__).with_name("results")
    names = (
        "quotient-reconstruction.json",
        "quotient-divisor-inputs.json",
        "quotient-divisor-verification.json",
        "symmetric-kummer-inversion-sage.json",
    )
    return {
        "status": "source_preflight_only_experiment_not_executed",
        "scope": {
            "calibration_matrix_fitted": False,
            "target_divisors_evaluated": False,
            "inverse_experiment_executed": False,
            "holdout_scalar_labels_computed": False,
        },
        "input_sha256": {
            name: sha256((results / name).read_bytes()).hexdigest() for name in names
        },
        "working_field_modulus_ascending": MODULUS,
        "source": {
            "p": 7,
            "extension_degree": 3,
            "modulus_ascending": [1, 0, 1, 1],
            "a": 248,
            "b": 42,
            "generator": [0, 176],
            "subgroup_order": 373,
            "kernel_enumeration_index": 153,
            "kernel_basis_rows": rows,
            "theta_seed": seed,
        },
        "pilot": {
            "known_scalar_labels": list(range(1, 33)),
            "full_section_vectors": 32,
            "full_section_columns": 64,
            "exact_sample_rank": 8,
            "normalized_full_section_matrix_sha256": digest.hexdigest(),
            "selected_section_indices": columns,
            "selected_section_rank_on_pilot": 8,
            "global_section_identity_proved_here": False,
        },
        "training_frame": {
            "known_scalar_labels": [i + 1 for i in frame],
            "point_coordinates": [points[i] for i in frame],
            "eight_coordinate_vectors_hex": [
                bytes(c for v in vectors[i] for c in coefficients(v)).hex()
                for i in frame
            ],
            "ninth_frame_point_coordinates_in_first_eight_basis_hex": bytes(
                c for v in weights for c in coefficients(v)
            ).hex(),
            "every_eight_of_nine_rank_eight": True,
            "target_frame_not_yet_checked": True,
        },
        "holdout_stream": {
            "tag": "ghs-calibration-holdout-v1",
            "counter_limit": 10000,
            "counters_consumed": counter,
            "challenge_count": len(challenges),
            "pilot_and_challenge_sign_classes_disjoint": True,
            "challenge_and_auxiliary_sign_classes_pairwise_disjoint": True,
        },
        "holdout_coordinate_inputs": challenges,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
