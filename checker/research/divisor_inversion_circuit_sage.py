"""Sparse ordered-support inversion attempt with explicit auxiliary variables.

This avoids expanding squared determinant/norm products in the symmetric
system. Four deterministic nonsingular, denominator-free support controls are
attempted; known supports validate equations but are not passed to the solver.
"""

import argparse
import itertools
import json
import multiprocessing
import queue
import sys
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from sage.all import GF, PolynomialRing, prod
from sage.env import SAGE_VERSION

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from research.divisor_inversion_system_sage import det3
from research.quotient_divisor_verification import normalize, verify_source

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def solver_worker(equations, A, result_queue):
    try:
        ideal = A.ideal(equations)
        basis = ideal.groebner_basis(algorithm="libsingular:std")
        dimension = int(ideal.dimension())
        result = {
            "status": "groebner_basis_computed",
            "groebner_basis_size": len(basis),
            "dimension": dimension,
        }
        if dimension == 0:
            solutions = ideal.variety()
            result["solution_count_in_working_field"] = len(solutions)
            result["solutions_first_six_coordinates"] = [
                [list(map(int, solution[A.gen(i)].polynomial())) for i in range(6)]
                for solution in solutions[:24]
            ]
        result_queue.put(result)
    except (
        RuntimeError,
        ArithmeticError,
        ValueError,
        TypeError,
        NotImplementedError,
    ) as error:
        result_queue.put(
            {
                "status": "solver_backend_error",
                "error": f"{type(error).__name__}: {error}",
            }
        )


def isolated_solver(equations, A, budget):
    context = multiprocessing.get_context("fork")
    result_queue = context.Queue()
    process = context.Process(target=solver_worker, args=(equations, A, result_queue))
    process.start()
    process.join(budget)
    if process.is_alive():
        process.terminate()
        process.join(5)
        if process.is_alive():
            process.kill()
            process.join()
        return {"status": "solver_time_budget_exhausted", "isolated_process": True}
    try:
        result = result_queue.get(timeout=1)
    except queue.Empty:
        result = {
            "status": "solver_process_exited_without_result",
            "exit_code": process.exitcode,
        }
    result["isolated_process"] = True
    return result


def build(ctx):
    k = ctx["field"]
    Kx = PolynomialRing(k, "t")
    basis = [ctx["phi"]] + ctx["psi"]
    denominators = []
    for functions in basis:
        den = Kx.one()
        for f in functions:
            for c in f.list():
                den = den.lcm(Kx(c.denominator()))
        denominators.append(den)
    norm_denominators = []
    for f in ctx["norm_functions"]:
        den = Kx.one()
        for c in f.list():
            den = den.lcm(Kx(c.denominator()))
        norm_denominators.append(den)
    names = [f"{v}{i}" for i in range(3) for v in ("x", "y")]
    names += [f"e{b}_{j}_{i}" for b in range(7) for j in range(3) for i in range(3)]
    names += [f"det{b}" for b in range(7)] + [f"den{b}" for b in range(7)]
    names += (
        [f"normn{i}" for i in range(6)]
        + [f"normd{i}" for i in range(6)]
        + ["saturation"]
    )
    A = PolynomialRing(k, names=names, order="degrevlex")
    vars = A.gens()
    xy = [vars[2 * i : 2 * i + 2] for i in range(3)]
    evaluations = [
        [[vars[6 + 9 * b + 3 * j + i] for i in range(3)] for j in range(3)]
        for b in range(7)
    ]
    dets, norms = vars[69:76], vars[76:83]
    normn, normd, saturation = vars[83:89], vars[89:95], vars[95]
    assert len(vars) == 96

    def eval_poly(poly, x):
        return sum(A(c) * x**i for i, c in enumerate(poly))

    def numerator(f, den, x, y):
        value = A.zero()
        for j, c in enumerate(f.list()):
            q, rem = den.quo_rem(Kx(c.denominator()))
            assert not rem
            value += eval_poly(q * Kx(c.numerator()), x) * y**j
        return value

    equations = []
    for x, y in xy:
        curve = ctx["curve"].defining_polynomial()
        equations.append(
            sum(A(c) * x ** mon[1] * y ** mon[2] for mon, c in curve.dict().items())
        )
    for b, functions in enumerate(basis):
        for j, f in enumerate(functions):
            for i, (x, y) in enumerate(xy):
                equations.append(
                    evaluations[b][j][i] - numerator(f, denominators[b], x, y)
                )
        equations.append(
            dets[b] - det3([[evaluations[b][j][i] for j in range(3)] for i in range(3)])
        )
        equations.append(norms[b] - prod(eval_poly(denominators[b], x) for x, y in xy))
    for j, f in enumerate(ctx["norm_functions"]):
        equations.append(
            normn[j] - prod(numerator(f, norm_denominators[j], x, y) for x, y in xy)
        )
        equations.append(
            normd[j] - prod(eval_poly(norm_denominators[j], x) for x, y in xy)
        )
    equations.append(saturation * dets[0] * prod(norms) * prod(normd) - 1)
    print(
        f"constructed {len(equations)} sparse equations in {len(vars)} variables",
        file=sys.stderr,
        flush=True,
    )

    def witness(support):
        out = [v for pair in support for v in pair]
        es, ds, ns, nns, nds = [], [], [], [], []
        for b, functions in enumerate(basis):
            entries = [
                [numerator(f, denominators[b], A(x), A(y)) for x, y in support]
                for f in functions
            ]
            entries = [[k(value) for value in row] for row in entries]
            es.extend(value for row in entries for value in row)
            ds.append(det3([[entries[j][i] for j in range(3)] for i in range(3)]))
            ns.append(prod(denominators[b](x) for x, y in support))
        for j, f in enumerate(ctx["norm_functions"]):
            nns.append(
                prod(
                    k(numerator(f, norm_denominators[j], A(x), A(y)))
                    for x, y in support
                )
            )
            nds.append(prod(norm_denominators[j](x) for x, y in support))
        product = ds[0] * prod(ns) * prod(nds)
        if not product:
            raise ValueError("zero chart denominator")
        out += es + ds + ns + nns + nds + [1 / product]
        assert len(out) == len(vars)
        target = [k(1)] + [
            ctx["beta"][j]
            * ds[j + 1] ** 2
            * ns[0] ** 2
            * nns[j]
            / (ds[0] ** 2 * ns[j + 1] ** 2 * nds[j])
            for j in range(6)
        ]
        return out, target

    def system(target):
        transfer = [
            ctx["beta"][j] * dets[j + 1] ** 2 * norms[0] ** 2 * normn[j]
            - target[j + 1] * dets[0] ** 2 * norms[j + 1] ** 2 * normd[j]
            for j in range(6)
        ]
        return equations + transfer

    return A, witness, system


def run(budget, output):
    root = Path(__file__).resolve().parents[1]
    path = root / "research/results/quotient-divisor-inputs.json"
    data = json.loads(path.read_text())
    start = perf_counter()
    ctx = verify_source(data, data["sources"][0], context_only=True)
    setup = perf_counter() - start
    before = perf_counter()
    A, witness, system = build(ctx)
    report = {
        "experiment": 4,
        "method": "ordered-support sparse auxiliary-variable system",
        "sage_version": SAGE_VERSION,
        "setup_seconds": setup,
        "construction_seconds": perf_counter() - before,
        "variables": A.ngens(),
        "input_sha256": sha256(path.read_bytes()).hexdigest(),
        "solver_budget_seconds_per_control": budget,
        "controls": [],
        "subgroup_enumeration_used": False,
        "scalable_inversion_established": False,
    }
    rational = []
    for P in ctx["base_curve"].rational_points():
        if P[0]:
            rational.append((ctx["field"](P[1] / P[0]), ctx["field"](P[2] / P[0])))
    L = GF(49, "w", modulus=PolynomialRing(GF(7), "t")([1, 0, 1]))
    w = ctx["field"](-1).sqrt()

    def embed(value):
        return sum(ctx["field"](c) * w**i for i, c in enumerate(value.polynomial()))

    for P in ctx["base_curve"].change_ring(L).rational_points():
        if P[0]:
            pair = (embed(P[1] / P[0]), embed(P[2] / P[0]))
            if pair not in rational:
                rational.append(pair)
    for support in itertools.combinations(rational, 3):
        if len({x for x, y in support}) != 3:
            continue
        try:
            values, target = witness(support)
        except ValueError:
            continue
        equations = system(target)
        assert all(eq(*values) == 0 for eq in equations)
        D = sum(ctx["curve"]([1, x, y]).place().divisor() for x, y in support)
        point = ctx["group"].point(D - ctx["origin_effective"])
        assert ctx["jacobian_image"](point) == normalize(target)
        control = {
            "support_for_validation_only_working_field_coefficients": [
                [list(map(int, x.polynomial())), list(map(int, y.polynomial()))]
                for x, y in support
            ],
            "equations": len(equations),
            "all_forward_equations_verified": True,
            "independent_Hess_class_evaluation_matches": True,
        }
        started = perf_counter()
        control.update(isolated_solver(equations, A, budget))
        control["solver_seconds"] = perf_counter() - started
        report["controls"].append(control)
        output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(
            f"sparse control {len(report['controls'])}: {control['status']}",
            file=sys.stderr,
            flush=True,
        )
        if len(report["controls"]) == 4:
            break
    assert len(report["controls"]) == 4
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget", type=int, default=120)
    args = parser.parse_args()
    result = run(args.budget, args.output)
    result["source_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
