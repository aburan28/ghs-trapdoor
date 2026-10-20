"""Measured toy BSGS baseline; backend differences preclude attack comparisons."""

import argparse
import json
from hashlib import sha256
from math import isqrt
from pathlib import Path
from time import perf_counter

from .small_fields import SmallField, ec_add, ec_mul

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def run():
    start = perf_counter()
    k = SmallField(7, 3, [1, 0, 1, 1])
    k.tables()
    a, generator, r = 248, (0, 176), 373
    inputs = Path(__file__).with_name("results") / "quotient-divisor-inputs.json"
    saved = json.loads(inputs.read_text())
    target = tuple(saved["sources"][0]["source_points_for_validation"][137])
    assert ec_mul(k, 0, a, generator, r) is None
    setup = perf_counter() - start
    before = perf_counter()
    m = isqrt(r) + 1
    table, point, additions = {}, None, 0
    for j in range(m):
        table[point] = j
        point = ec_add(k, 0, a, point, generator)
        additions += 1
    step = (point[0], k.neg(point[1]))
    giant = target
    solution = None
    for j in range(m + 1):
        if giant in table:
            solution = (j * m + table[giant]) % r
            break
        giant = ec_add(k, 0, a, giant, step)
        additions += 1
    evaluation = perf_counter() - before
    assert solution == 137 and ec_mul(k, 0, a, generator, solution) == target
    return {
        "field": "F343",
        "subgroup_order": r,
        "baseline": "baby-step giant-step on the original elliptic subgroup",
        "input_generator": generator,
        "input_target": target,
        "solution_for_validation": solution,
        "baby_table_entries": len(table),
        "counted_group_additions": additions,
        "timing_seconds": {"setup": setup, "DLP_evaluation": evaluation},
        "input_sha256": sha256(inputs.read_bytes()).hexdigest(),
        "unimplemented_attack_cost_categories": {
            "trace_zero_relation_generation": None,
            "linear_algebra": None,
            "large_field_divisor_inversion": None,
        },
        "comparison_scope": "toy sanity baseline only; different arithmetic backends and representations; no end-to-end speedup comparison",
        "practical_speedup_established": False,
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
