"""Exhaust every level-two Lagrangian in declared dimension-four modules.

No Jacobian-locus conclusion follows from this necessary kernel screen.
"""

import argparse
import itertools
import json
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from ghs_checker.analysis import binary_basis

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")

def rref_four_planes():
    for pivots in itertools.combinations(range(8), 4):
        free = [
            (i, j) for i in range(4) for j in range(pivots[i] + 1, 8) if j not in pivots
        ]
        for mask in range(1 << len(free)):
            rows = [1 << p for p in pivots]
            for bit, (i, j) in enumerate(free):
                if mask >> bit & 1:
                    rows[i] |= 1 << j
            yield tuple(rows)


def transform(v, rows):
    out = 0
    for i, row in enumerate(rows):
        if v >> i & 1:
            out ^= row
    return out


def contains(v, rows):
    for row in rows:
        if v & (row & -row):
            v ^= row
    return not v


def elliptic_form(v):
    return ((v & 0b01010101) << 1) | ((v & 0b10101010) >> 1)


def genus_two_form(v):
    out = 0
    for shift in (0, 4):
        block = (v >> shift) & 15
        # Even branch subsets modulo complement, vi={i,5}: <vi,vj>=1 for i!=j.
        parity = block.bit_count() & 1
        out |= (block ^ (15 if parity else 0)) << shift
    return out


def partitions(n, minimum=1):
    if n == 0:
        yield ()
    for first in range(minimum, n+1):
        for rest in partitions(n-first, first):
            yield (first,)+rest


def branch_permutation_module(partition):
    permutation = list(range(6))
    offset = 0
    for length in partition:
        for j in range(length):
            permutation[offset+j] = offset+(j+1) % length
        offset += length
    rows = []
    for i in range(4):
        subset = (1 << permutation[i]) ^ (1 << permutation[5])
        if subset & 16:
            subset ^= 63  # quotient by the full six-point branch subset
        rows.append(subset & 15)
    return rows


def run():
    start = perf_counter()
    elliptic, genus_two, planes = [], [], 0
    for rows in rref_four_planes():
        planes += 1
        for target, form in ((elliptic, elliptic_form), (genus_two, genus_two_form)):
            if all(
                (form(a) & b).bit_count() % 2 == 0
                for a, b in itertools.combinations(rows, 2)
            ):
                target.append(rows)
    assert planes == 200787 and len(elliptic) == len(genus_two) == 2295
    setup = perf_counter() - start
    modules = []
    matrices = {
        "identity": (1, 2),
        "unipotent": (2, 1),
        "irreducible_order3": (2, 3),
    }
    for name, block in matrices.items():
        frob = [1 << (j + 2) for j in range(6)] + list(block)
        modules.append(
            (
                f"elliptic_degree4_{name}",
                elliptic,
                elliptic_form,
                frob,
                [(2 * i, 2) for i in range(4)],
            )
        )
    cycle_types = {}
    for partition in partitions(6):
        name = "genus2_degree2_split_2torsion" if partition == (1,)*6 else "genus2_degree2_cycles_"+"_".join(map(str, partition))
        cycle_types[name] = partition
        modules.append((name, genus_two, genus_two_form,
                        [1 << (j+4) for j in range(4)]+branch_permutation_module(partition),
                        [(0, 4), (4, 4)]))
    assert len(cycle_types) == 11
    output = []
    for name, candidates, form, frob, blocks in modules:
        # Verify Frobenius preserves the alternating form on all basis pairs.
        for i, j in itertools.product(range(8), repeat=2):
            assert (form(1 << i) & (1 << j)).bit_count() % 2 == (
                form(frob[i]) & frob[j]
            ).bit_count() % 2
        stable = []
        for rows in candidates:
            if not all(contains(transform(row, frob), rows) for row in rows):
                continue
            intersections = []
            for shift, dimension in blocks:
                mask = ((1 << dimension) - 1) << shift
                intersections.append(
                    4 - len(binary_basis([row & (255 ^ mask) for row in rows]))
                )
            # Independent all-vector check, not only basis stability.
            for subset in range(16):
                v = 0
                for i in range(4):
                    if subset >> i & 1:
                        v ^= rows[i]
                assert contains(transform(v, frob), rows)
            stable.append(
                {
                    "basis_rows_binary": rows,
                    "factor_intersection_dimensions": intersections,
                    "passes_no_factor_test": not any(intersections),
                }
            )
        output.append(
            {
                "module": name,
                "Frobenius_rows_binary": frob,
                "all_geometric_Lagrangians": len(candidates),
                "rational_Lagrangians": len(stable),
                "passing_no_factor_test": sum(
                    row["passes_no_factor_test"] for row in stable
                ),
                "kernels": stable,
                "branch_cycle_partition": cycle_types.get(name),
            }
        )
    return {
        "experiment": 6,
        "ell": 2,
        "dimension": 4,
        "RREF_planes_examined_per_pairing": planes,
        "modules": output,
        "timing_seconds": {
            "setup": setup,
            "evaluation": perf_counter() - start - setup,
        },
        "scope": {
            "all_declared_level2_module_kernels_exhausted": True,
            "module_curve_realization_verified_in_this_script": False,
            "principally_polarized_quotient_kernel_condition_verified": True,
            "indecomposability_certified": False,
            "Jacobian_locus_certified": False,
            "Schottky_filter_implemented": False,
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
