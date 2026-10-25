"""Compare normalized two-pole class genera with the existing GHS checker.

The class sweep is over a fixed x-coordinate. Scaling x identifies many
presentations. Geometric class-span genus does not remove the checker's
base-field regularity condition.
"""

import argparse
import json
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from ghs_checker.analysis import analyze_descent
from ghs_checker.curve import BinaryCurve
from ghs_checker.field import BinaryField

from .wild_artin_schreier import class_span, span

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def run():
    fields = []
    for n, modulus in ((4, 19), (6, 67)):
        start = perf_counter()
        k = BinaryField(n, modulus)
        q = k.size
        trace = [k.absolute_trace(x) for x in range(q)]
        constant = next(x for x in range(q) if trace[x])
        rows, regular, irregular, maps = [], 0, 0, 0
        for beta in range(1, q):
            for bit in (0, 1):
                c = constant if bit else 0
                E = BinaryCurve(k, c, k.square(beta))
                # y_AS^2+y_AS=x+beta/x+c;
                # (X,Y)=(x,x*y_AS+beta) is the standard binary model.
                count = 2
                for x in range(1, q):
                    rhs = x ^ k.div(beta, x) ^ c
                    for y in range(q):
                        if k.square(y) ^ y == rhs:
                            assert E.contains((x, k.mul(x, y) ^ beta))
                            count += 1
                            maps += 1
                assert E.count_points() == count
                for d in range(1, n):
                    if n % d:
                        continue
                    basis = class_span(k, 1, beta, bit, subfield_degree=d)
                    values = span(basis)
                    assert (1 << (2 * n)) not in values
                    genus = sum(
                        bool(v & (q - 1)) and bool((v >> n) & (q - 1)) for v in values
                    )
                    baseline = analyze_descent(E, d)
                    assert len(basis) == baseline["m"]
                    if baseline["regularity_condition_holds"]:
                        assert genus == baseline["genus"]
                        regular += 1
                    else:
                        irregular += 1
                    rows.append(
                        {
                            "beta": beta,
                            "constant_trace_bit": bit,
                            "subfield_degree": d,
                            "source_order": count,
                            "class_span_rank": len(basis),
                            "geometric_span_genus_over_K": genus,
                            "classical_regularity": baseline[
                                "regularity_condition_holds"
                            ],
                            "classical_GHS_genus": baseline["genus"],
                            "classical_epsilon": baseline["one_zero_in_augmented_span"],
                        }
                    )
        fields.append(
            {
                "degree": n,
                "modulus_binary": modulus,
                "raw_nonzero_two_pole_presentations": 2 * (q - 1) ** 2,
                "x_scaling_normalized_presentations": 2 * (q - 1),
                "regular_comparisons_equal_genus": regular,
                "irregular_cases_not_claimed_as_base_field_descents": irregular,
                "affine_model_point_maps_checked": maps,
                "comparisons": rows,
                "elapsed_seconds": perf_counter() - start,
            }
        )
    return {
        "experiment": 8,
        "fields": fields,
        "normalization": "X=a*x transforms (a,b,c) to (1,a*b,c); then Y=X*y_AS+a*b gives BinaryCurve(a=c,b=(a*b)^2)",
        "scope": {
            "all_declared_normalized_ordinary_cases_compared": True,
            "lower_genus_than_regular_classical_baseline_found": False,
            "irregularity_condition_bypassed": False,
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
