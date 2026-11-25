"""Exhaustive reduced one/two-pole Artin--Schreier class experiments.

Run from checker/: python -m research.wild_artin_schreier --output FILE
The span genus is a conductor/discriminant calculation, not a constructed
descent curve or an evaluated Jacobian subgroup map.
"""

import argparse
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from ghs_checker.analysis import binary_basis
from ghs_checker.field import BinaryField

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")

def span(basis):
    values = [0]
    for b in basis:
        values += [v ^ b for v in values]
    return values


def class_span(k, a, b, constant_bit, subfield_degree=1):
    conjugates = []
    for _ in range(k.degree//subfield_degree):
        conjugates.append(a | (b << k.degree) | (constant_bit << (2 * k.degree)))
        a, b = k.frobenius(a, subfield_degree), k.frobenius(b, subfield_degree)
    basis = binary_basis(conjugates)
    return [basis[pivot] for pivot in sorted(basis, reverse=True)]


def run_field(n, modulus):
    start = perf_counter()
    k = BinaryField(n, modulus)
    q = k.size
    trace = [k.absolute_trace(x) for x in range(q)]
    mul = [[k.mul(a, x) for x in range(q)] for a in range(q)]
    cubes = [k.pow(x, 3) for x in range(q)]
    inverses = [0] + [k.inv(x) for x in range(1, q)]
    constant = next(x for x in range(q) if trace[x])
    image_wp = {k.square(h) ^ h for h in range(q)}
    assert image_wp == {x for x in range(q) if trace[x] == 0}
    assert len(image_wp) == q // 2
    setup = perf_counter() - start
    families = []
    for family in ("cubic", "two_poles"):
        digest, stats, examples = sha256(), Counter(), {}
        inspected = geometric = 0
        for a in range(q):
            for b in range(q):
                for bit in (0, 1):
                    inspected += 1
                    basis = class_span(k, a, b, bit)
                    values = span(basis)
                    assert len(values) == 2 ** len(basis)
                    # A pure nontrivial constant creates a constant extension.
                    connected = (1 << (2 * n)) not in values
                    aa = [v & (q - 1) for v in values]
                    bb = [(v >> n) & (q - 1) for v in values]
                    genus = (
                        sum(x != 0 for x in aa)
                        if family == "cubic"
                        else sum(x != 0 and y != 0 for x, y in zip(aa, bb))
                    )
                    if family == "cubic":
                        count = (
                            1
                            + 2
                            * sum(
                                trace[mul[a][cubes[x]] ^ mul[b][x]] == bit
                                for x in range(q)
                            )
                            if a or b
                            else None
                        )
                        source_genus = int(a != 0)
                    else:
                        domain = range(1, q) if b else range(q)
                        poles = int(a != 0) + int(b != 0)
                        count = (
                            poles
                            + 2
                            * sum(
                                trace[mul[a][x] ^ mul[b][inverses[x]]] == bit
                                for x in domain
                            )
                            if a or b
                            else None
                        )
                        source_genus = int(a != 0 and b != 0)
                    if source_genus:
                        t = q + 1 - count
                        assert t * t <= 4 * q
                        polynomial = [q, -t, 1]
                    else:
                        polynomial = [1] if a or b else None
                    stats[(len(basis), connected, genus, source_genus)] += 1
                    geometric += connected and bool(a or b)
                    row = [a, b, bit, count, len(basis), connected, genus, polynomial]
                    digest.update(
                        (json.dumps(row, separators=(",", ":")) + "\n").encode()
                    )
                    key = (len(basis), connected, genus, source_genus)
                    if key not in examples and a and b:
                        examples[key] = {
                            "a": a,
                            "b": b,
                            "constant_trace_bit": bit,
                            "constant_representative": constant if bit else 0,
                            "Frobenius_span_basis_binary": basis,
                            "span_rank": len(basis),
                            "geometrically_connected": connected,
                            "compositum_genus_from_conductors": genus
                            if connected
                            else None,
                            "source_genus": source_genus,
                            "source_point_count": count,
                            "source_full_Frobenius_polynomial_ascending": polynomial,
                        }
        families.append(
            {
                "family": family,
                "reduced_classes_examined": inspected,
                "geometric_nonconstant_spans": geometric,
                "ordered_class_results_sha256": digest.hexdigest(),
                "distribution": [
                    {
                        "rank": r,
                        "geometric": c,
                        "genus": g if c else None,
                        "source_genus": sg,
                        "classes": number,
                    }
                    for (r, c, g, sg), number in sorted(stats.items())
                ],
                "representatives": list(examples.values()),
                "ordinary_genus1_family": family == "two_poles",
            }
        )
    # All even-pole coefficients reduce by adjoining h=sqrt(c)*x or /x.
    for c in range(q):
        root = k.pow(c, q // 2)
        assert k.square(root) == c
        for x in range(q):
            h = mul[root][x]
            assert mul[c][k.square(x)] ^ mul[root][x] == k.square(h) ^ h
    return {
        "degree": n,
        "modulus_binary": modulus,
        "field_size": q,
        "constant_wp_image_size": len(image_wp),
        "constant_trace_one_representative": constant,
        "families": families,
        "timing_seconds": {
            "setup": setup,
            "evaluation": perf_counter() - start - setup,
        },
    }


def run():
    return {
        "experiment": 8,
        "bounds": {
            "fields": [[4, 0b10011], [6, 0b1000011]],
            "coefficients": "all a,b and both constant trace classes",
            "families": ["a*x^3+b*x+c", "a*x+b/x+c"],
        },
        "fields": [run_field(4, 0b10011), run_field(6, 0b1000011)],
        "scope": {
            "all_declared_reduced_classes_exhausted": True,
            "full_Frobenius_polynomials_for_genus1_quotients": True,
            "compositum_genus_uses_conductor_discriminant_formula": True,
            "general_compositum_equations_constructed": False,
            "Jacobian_norm_conorm_evaluated": False,
            "practical_speedup_established": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run()
    report["source_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
