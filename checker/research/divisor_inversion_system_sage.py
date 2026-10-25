"""Construct and attempt the symmetric degree-three Weil inversion system.

There is no subgroup table in this solver. Bounds are solver wall-clock
budgets, not mathematical nonexistence claims. Run with sage -python.
"""

import argparse
import json
import sys
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from cysignals.alarm import AlarmInterrupt, alarm, cancel_alarm
from sage.all import PolynomialRing, prod
from sage.env import SAGE_VERSION

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from research.quotient_divisor_verification import verify_source

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")

def remainder(f, u):
    x = f.parent().gen()
    while f and f.degree() >= 3:
        f -= f.leading_coefficient() * x ** (f.degree() - 3) * u
    return f


def det3(a):
    # Avoid fraction-field elimination and its multivariate gcds.
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def construct(ctx):
    k = ctx["field"]
    A = PolynomialRing(k, names=("u0", "u1", "u2", "v0", "v1", "v2"), order="degrevlex")
    u0, u1, u2, v0, v1, v2 = A.gens()
    R = PolynomialRing(A, "x")
    x = R.gen()
    U, V = x**3 + u2 * x * x + u1 * x + u0, v2 * x * x + v1 * x + v0
    Kx = PolynomialRing(k, "t")

    def convert(f):
        return R([A(c) for c in f])

    def polynomial_function(f, denominator):
        result = R.zero()
        for j, c in enumerate(f.list()):
            numerator, den = Kx(c.numerator()), Kx(c.denominator())
            quotient, rem = denominator.quo_rem(den)
            assert rem == 0
            result += convert(numerator * quotient) * V**j
        return remainder(result, U)

    def common(functions):
        denominator = Kx.one()
        for f in functions:
            for c in f.list():
                denominator = denominator.lcm(Kx(c.denominator()))
        return denominator

    def norm(f):
        columns = [remainder(x**j * f, U) for j in range(3)]
        return det3([[columns[j][i] for j in range(3)] for i in range(3)])

    def determinant(functions):
        denominator = common(functions)
        columns = [polynomial_function(f, denominator) for f in functions]
        value = det3([[columns[j][i] for j in range(3)] for i in range(3)])
        return value, norm(convert(denominator))

    start = perf_counter()
    dphi, nphi = determinant(ctx["phi"])
    print("phi determinant constructed", file=sys.stderr, flush=True)
    numerators, denominators = [], []
    for j, (functions, f) in enumerate(zip(ctx["psi"], ctx["norm_functions"])):
        dpsi, npsi = determinant(functions)
        print(f"psi determinant {j} constructed", file=sys.stderr, flush=True)
        den = common([f])
        nf = norm(polynomial_function(f, den))
        df = norm(convert(den))
        numerator = ctx["beta"][j] * dpsi**2 * nphi**2 * nf
        denominator = dphi**2 * npsi**2 * df
        # Cancel polynomial artifacts before adding the target values.
        if max(numerator.total_degree(), denominator.total_degree()) <= 20:
            divisor = numerator.gcd(denominator)
            numerator, denominator = numerator // divisor, denominator // divisor
        numerators.append(numerator)
        denominators.append(denominator)
        print(
            f"Weil equation {j}: degrees {numerator.total_degree()}/{denominator.total_degree()}",
            file=sys.stderr,
            flush=True,
        )
    curve = ctx["curve"].defining_polynomial()
    equation = sum(A(c) * x ** mon[1] * V ** mon[2] for mon, c in curve.dict().items())
    curve_remainder = remainder(equation, U)
    return (
        A,
        numerators,
        denominators,
        [curve_remainder[j] for j in range(3)],
        perf_counter() - start,
    )


def chart_coordinates(ctx, point):
    D = point.defining_divisor()
    if any(m != 1 for p, m in D.list()):
        raise ValueError("repeated support outside this chart")
    _, high, lift = ctx["evaluations"](D)
    L, _, _ = high.residue_field()
    pairs = []
    for place, _ in D.list():
        xx, yy = [f.evaluate(place) for f in ctx["coordinates"][1:]]
        if place.degree() == 1 and high.degree() > 1:
            xx, yy = lift(xx), lift(yy)
        for j in range(place.degree()):
            pairs.append(
                (
                    xx ** (ctx["field"].cardinality() ** j),
                    yy ** (ctx["field"].cardinality() ** j),
                )
            )
    if len({p[0] for p in pairs}) != 3:
        raise ValueError("projection does not separate support")
    R = PolynomialRing(L, "x")
    x = R.gen()
    U = prod(x - a for a, b in pairs)
    V = R.lagrange_polynomial(pairs)
    return ctx["descend"]([U[j] for j in range(3)] + [V[j] for j in range(3)], high)


def run(budget):
    root = Path(__file__).resolve().parents[1]
    path = root / "research/results/quotient-divisor-inputs.json"
    data = json.loads(path.read_text())
    start = perf_counter()
    ctx = verify_source(data, data["sources"][0], context_only=True)
    setup = perf_counter() - start
    result = {
        "experiment": 4,
        "sage_version": SAGE_VERSION,
        "input_sha256": sha256(path.read_bytes()).hexdigest(),
        "setup_seconds": setup,
        "solver_budget_seconds_per_attempt": budget,
        "chart": "U=x^3+u2*x^2+u1*x+u0, Z=V(x), V=v2*x^2+v1*x+v0",
        "subgroup_enumeration_used": False,
        "scope": {
            "scalable_inversion_established": False,
            "global_sign_resolved": False,
        },
        "attempts": [],
    }
    # Construction also has an explicit budget; a timeout is an engineering result.
    try:
        alarm(budget)
        A, numerators, denominators, curve_equations, construction = construct(ctx)
        cancel_alarm()
    except AlarmInterrupt:
        cancel_alarm()
        result["construction_status"] = "time_budget_exhausted"
        return result
    result["construction_status"] = "complete"
    result["construction_seconds"] = construction
    result["equations"] = [
        {
            "numerator_degree": int(n.total_degree()),
            "denominator_degree": int(d.total_degree()),
            "numerator_terms": len(n.monomials()),
            "denominator_terms": len(d.monomials()),
            "numerator_sha256": sha256(str(n).encode()).hexdigest(),
            "denominator_sha256": sha256(str(d).encode()).hexdigest(),
        }
        for n, d in zip(numerators, denominators)
    ]
    for scalar in (1, 2, 37, 111):
        attempt = {"validation_scalar_not_given_to_solver": scalar}
        point = scalar * ctx["generator"]
        target = ctx["jacobian_image"](point)
        try:
            coordinates = chart_coordinates(ctx, point)
        except ValueError as error:
            attempt["status"] = str(error)
            result["attempts"].append(attempt)
            continue
        if not target[0] or any(d(*coordinates) == 0 for d in denominators):
            attempt["status"] = "unsupported_zero_denominator_chart"
            result["attempts"].append(attempt)
            continue
        assert all(f(*coordinates) == 0 for f in curve_equations)
        equations = [
            n - (value / target[0]) * d
            for n, d, value in zip(numerators, denominators, target[1:])
        ]
        assert all(f(*coordinates) == 0 for f in equations)
        attempt["forward_polynomial_identity_verified"] = True
        B = PolynomialRing(
            ctx["field"], names=A.variable_names() + ("saturation",), order="degrevlex"
        )
        saturation = B.gen(6) * B(prod(denominators)) - 1
        ideal = B.ideal([B(f) for f in curve_equations + equations] + [saturation])
        before = perf_counter()
        try:
            alarm(budget)
            basis = ideal.groebner_basis()
            attempt["groebner_basis_length"] = len(basis)
            attempt["dimension"] = int(ideal.dimension())
            # Do not feed known coordinates to variety(), even for validation.
            solutions = ideal.variety() if ideal.dimension() == 0 else []
            attempt["solutions_in_working_field"] = len(solutions)
            attempt["status"] = (
                "finite_solutions_obtained"
                if solutions
                else "no_working_field_solutions_returned"
            )
            cancel_alarm()
        except AlarmInterrupt:
            cancel_alarm()
            attempt["status"] = "groebner_time_budget_exhausted"
        attempt["solver_seconds"] = perf_counter() - before
        result["attempts"].append(attempt)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--budget", type=int, default=120)
    args = parser.parse_args()
    result = run(args.budget)
    result["source_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
