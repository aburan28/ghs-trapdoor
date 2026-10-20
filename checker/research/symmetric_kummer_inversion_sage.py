"""Sage: nonenumerative inversion of a genus-3 symmetric Kummer model.

Construct the 8-dimensional diagonal-vanishing system in Tian section 3.2,
then invert supplied projective values by a zero-dimensional polynomial
system. The linear calibration from the elliptic quotient's theta sections
to this basis is NOT implemented here. No Jacobian group is enumerated.
"""

import argparse
import itertools
import json
from hashlib import sha256
from pathlib import Path

from sage.all import GF, Curve, PolynomialRing, matrix, prod, vector
from sage.env import SAGE_VERSION

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def coefficients(value, degree):
    return [int(value.polynomial()[i]) for i in range(degree)]


def basis_functions(x, y):
    return [1, x, y, x * x, x * y, y * y, x**3, x * x * y]


def triangular_solutions(polynomials, field, stats, known=None):
    """Bounded-degree roots and constant-coefficient linear elimination.

    This handles the triangular steps occurring in these test fibers; it is
    not a general replacement for a computer algebra system's variety solver.
    """
    known = {} if known is None else known
    polynomials = [g for g in polynomials if g != 0]
    if any(g.is_constant() for g in polynomials):
        return []
    if not polynomials:
        return [known]
    ring = polynomials[0].parent()
    for polynomial in polynomials:
        for variable in polynomial.variables():
            if polynomial.degree(variable) != 1:
                continue
            coefficient = polynomial.derivative(variable)
            if not coefficient.is_constant():
                continue
            expression = -(polynomial - coefficient * variable) / coefficient
            substitutions = {variable: expression}
            updated = {v: g.subs(substitutions) for v, g in known.items()}
            updated[variable] = expression
            stats["linear_eliminations"] += 1
            return triangular_solutions(
                [g.subs(substitutions) for g in polynomials], field, stats, updated
            )
    univariate = sorted(
        (g for g in polynomials if len(g.variables()) == 1),
        key=lambda g: g.total_degree(),
    )
    if not univariate:
        raise ArithmeticError(
            "Fiber needs a triangular step not implemented in this bounded solver"
        )
    polynomial = univariate[0]
    variable = polynomial.variables()[0]
    uring = PolynomialRing(field, "t")
    unipoly = uring(
        [
            polynomial.monomial_coefficient(variable**i)
            for i in range(polynomial.degree(variable) + 1)
        ]
    )
    stats["univariate_factorizations"] += 1
    stats["maximum_univariate_degree"] = max(
        stats["maximum_univariate_degree"], int(unipoly.degree())
    )
    solutions = []
    for root in unipoly.roots(multiplicities=False):
        substitutions = {variable: ring(root)}
        updated = {v: g.subs(substitutions) for v, g in known.items()}
        updated[variable] = ring(root)
        solutions.extend(
            triangular_solutions(
                [g.subs(substitutions) for g in polynomials], field, stats, updated
            )
        )
    return solutions


def transformed_quartic(original, plane):
    """Find the rational bitangent and put its two contacts at infinity."""
    base = GF(7)
    field = plane.base_ring()
    X, Y, Z = plane.gens()
    uring = PolynomialRing(field, "t")
    t = uring.gen()
    lines = [(1, a, b) for a, b in itertools.product(base, repeat=2)]
    lines += [(0, 1, b) for b in base] + [(0, 0, 1)]
    for line in lines:
        pivot = next(i for i, c in enumerate(line) if c)
        free = [i for i in range(3) if i != pivot]
        forms = [uring(0)] * 3
        forms[free[0]], forms[free[1]] = t, uring(1)
        forms[pivot] = -sum(line[i] * forms[i] for i in free) / line[pivot]
        binary = sum(
            c * prod(forms[i] ** m[i] for i in range(3))
            for m, c in original.dict().items()
        )
        if any(n % 2 for _, n in binary.factor()) or (4 - binary.degree()) % 2:
            continue
        contacts = [
            vector(field, [h(root) for h in forms])
            for root in binary.roots(multiplicities=False)
        ]
        if binary.degree() < 4:
            point = [field(0)] * 3
            point[free[0]], point[pivot] = field(1), -line[free[0]] / line[pivot]
            contacts.append(vector(field, point))
        if len(contacts) != 2:
            continue
        P, Q = contacts
        T = vector(field, [int(i == pivot) for i in range(3)])
        change = matrix([P, Q, T]).transpose()
        assert change.det() != 0
        moved = original(*list(change * vector(plane, (X, Y, Z))))
        assert moved(X, Y, 0) == moved.monomial_coefficient(X**2 * Y**2) * X**2 * Y**2
        return moved, line, change
    raise ArithmeticError("No supported rational bitangent with two split contacts")


def run(path):
    data = json.loads(path.read_text())
    row = data["sources"][0]["targets"][0]
    Fp = GF(7)
    univariate = PolynomialRing(Fp, "t")
    field = GF(49, "w", modulus=univariate([1, 0, 1]))
    plane = PolynomialRing(field, names=("X", "Y", "Z"))
    X, Y, Z = plane.gens()
    original = sum(
        field(c) * X ** m[0] * Y ** m[1] * Z ** m[2]
        for m, c in row["quartic_coefficients_ascending_monomials"]
    )
    quartic, bitangent, change = transformed_quartic(original, plane)
    curve = Curve(quartic)
    affine = sorted(
        (tuple(P[:2]) for P in curve.rational_points() if P[2] == 1),
        key=lambda P: tuple(c for v in P for c in coefficients(v, 2)),
    )
    triples = list(itertools.combinations_with_replacement(range(8), 3))
    values = [basis_functions(x, y) for x, y in affine]
    rows = []
    for a in values[:32]:
        for b in values[:24]:
            rows.append(
                [
                    2 * (a[i] * a[j] * b[k] + a[i] * a[k] * b[j] + a[j] * a[k] * b[i])
                    for i, j, k in triples
                ]
            )
    diagonal = matrix(field, rows)
    kernel = diagonal.right_kernel().basis_matrix()
    assert diagonal.rank() == 112 and kernel.nrows() == 8
    assert all(c**7 == c for row in kernel for c in row)
    ring = PolynomialRing(
        field, names=("x1", "y1", "x2", "y2", "x3", "y3", "h"), order="degrevlex"
    )
    x1, y1, x2, y2, x3, y3, h = ring.gens()
    basis = [basis_functions(x, y) for x, y in ((x1, y1), (x2, y2), (x3, y3))]
    symmetric = [
        sum(
            prod(basis[i][indices[i]] for i in range(3))
            for indices in itertools.permutations(ijk)
        )
        for ijk in triples
    ]
    zeta = [sum(c * f for c, f in zip(row, symmetric)) for row in kernel]
    # Prove the diagonal identities as polynomial identities on C x C,
    # rather than relying only on a finite evaluation sample.
    diagonal_ideal = ring.ideal([quartic(x1, y1, 1), quartic(x3, y3, 1)])
    for function in zeta:
        assert diagonal_ideal.reduce(function.subs({x2: x1, y2: y1})) == 0

    # A fixed extension of degree six over F49 contains every degree<=3
    # residue field of either recovered effective divisor.
    splitting = GF(7**12, "v")
    sr = PolynomialRing(splitting, "t")
    image_w = min(
        sr([1, 0, 1]).roots(multiplicities=False),
        key=lambda a: tuple(coefficients(a, 12)),
    )
    embedding = field.hom([image_w], splitting)
    wide = ring.change_ring(splitting)
    wide_gens = wide.gens()

    def lift(polynomial):
        return wide({m: embedding(c) for m, c in polynomial.dict().items()})

    cases, used = [], set()
    counter = 0
    while len(cases) < 3:
        digest = sha256(f"ghs-symmetric-inversion-v1:{counter}".encode()).digest()
        counter += 1
        indices = sorted({digest[i] % len(affine) for i in range(3)})
        if len(indices) != 3 or tuple(indices) in used:
            continue
        points = [affine[i] for i in indices]
        if (
            len({x for x, _ in points}) != 3
            or matrix([[x, y, 1] for x, y in points]).det() == 0
        ):
            continue
        used.add(tuple(indices))
        arguments = [v for P in points for v in P] + [field(1)]
        target = [f(*arguments) for f in zeta]
        pivot = next(i for i, c in enumerate(target) if c != 0)
        equations = [ring(quartic(x, y, 1)) for x, y in ((x1, y1), (x2, y2), (x3, y3))]
        equations += [
            target[pivot] * zeta[i] - target[i] * zeta[pivot]
            for i in range(8)
            if i != pivot
        ]
        equations += [h * zeta[pivot] - 1]
        ideal = ring.ideal(equations)
        groebner = ideal.groebner_basis()
        degree = len(ideal.normal_basis())
        assert ideal.dimension() == 0 and degree == 12
        stats = {
            "linear_eliminations": 0,
            "univariate_factorizations": 0,
            "maximum_univariate_degree": 0,
        }
        solutions = triangular_solutions([lift(g) for g in groebner], splitting, stats)
        assert len(solutions) == 12
        unordered = {}
        for solution in solutions:
            coordinates = [
                splitting(solution[v].constant_coefficient()) for v in wide_gens
            ]
            assert all(lift(f)(*coordinates) == 0 for f in equations)
            support = tuple(
                sorted(
                    (tuple(coordinates[2 * i : 2 * i + 2]) for i in range(3)),
                    key=lambda P: tuple(c for v in P for c in coefficients(v, 12)),
                )
            )
            unordered[support] = unordered.get(support, 0) + 1
        assert sorted(unordered.values()) == [6, 6]
        original_support = tuple(
            sorted(
                (tuple(embedding(c) for c in P) for P in points),
                key=lambda P: tuple(c for v in P for c in coefficients(v, 12)),
            )
        )
        assert original_support in unordered
        divisors = list(unordered)
        # Exact residual-divisor certificate: D+Dbar is cut by a section of
        # L(K_C+2*P_infinity), whose basis here is {1,x,y,x^2}.
        conic_matrix = matrix(
            splitting, [[1, x, y, x * x] for D in divisors for x, y in D]
        )
        conics = conic_matrix.right_kernel().basis()
        assert len(conics) == 1
        conic = conics[0]
        assert conic[2] != 0
        tt = sr.gen()
        yy = -(conic[0] + conic[1] * tt + conic[3] * tt * tt) / conic[2]
        residual = sr(
            sum(
                embedding(c) * tt ** m[0] * yy ** m[1]
                for m, c in quartic.dict().items()
            )
        )
        expected = prod(tt - x for D in divisors for x, _ in D)
        assert (
            residual.degree() == 6
            and residual / residual.leading_coefficient() == expected
        )
        cases.append(
            {
                "input_divisor_points_F49": [
                    [coefficients(v, 2) for v in P] for P in points
                ],
                "input_eight_projective_values_F49": [
                    coefficients(v, 2) for v in target
                ],
                "zero_dimensional_fiber_degree": degree,
                "ordered_solutions": len(solutions),
                "unordered_effective_divisors": len(unordered),
                "input_divisor_recovered": True,
                "all_equations_verified_at_all_solutions": True,
                "six_permutations_per_divisor": True,
                "residual_relation": "D+Dbar ~ K_C+2*P_infinity",
                "residual_relation_certified_by_degree_6_polynomial_identity": True,
                "residual_conic_coefficients_F7_12": [
                    coefficients(c, 12) for c in conic
                ],
                "recovered_divisor_supports_F7_12": [
                    [[coefficients(v, 12) for v in P] for P in D] for D in divisors
                ],
                "solver_steps": stats,
            }
        )
    assert all(
        f(*([affine[0][0], affine[0][1]] * 2 + list(affine[1]) + [1])) == 0
        for f in zeta
    )
    return {
        "scope": {
            "toy_only": True,
            "Sage_version_used": SAGE_VERSION,
            "Jacobian_group_enumeration_used": False,
            "scalar_labels_used_for_inversion": False,
            "input_is_a_projective_point_in_the_constructed_symmetric_Kummer_model": True,
            "calibration_from_elliptic_quotient_theta_coordinates_implemented": False,
            "end_to_end_nonenumerative_EC_to_Jacobian_transfer": False,
            "exceptional_divisors_with_repeated_support_or_infinity_handled": False,
            "general_triangular_solver_implemented": False,
            "practical_DLP_speedup_established": False,
        },
        "input_reconstruction_sha256": sha256(path.read_bytes()).hexdigest(),
        "source_kernel_index": row["kernel_enumeration_index"],
        "bitangent_over_F7": list(map(int, bitangent)),
        "change_of_plane_coordinates_F49": [
            [coefficients(c, 2) for c in row] for row in change
        ],
        "transformed_quartic": [
            [list(m), coefficients(c, 2)] for m, c in sorted(quartic.dict().items())
        ],
        "basis_functions": ["1", "x", "y", "x^2", "x*y", "y^2", "x^3", "x^2*y"],
        "symmetric_products": 120,
        "symmetrization_includes_all_six_permutations": True,
        "diagonal_sample_matrix_rows": diagonal.nrows(),
        "diagonal_matrix_rank": 112,
        "diagonal_vanishing_space_dimension": 8,
        "diagonal_identities_verified_symbolically": True,
        "zeta_total_degrees": [int(f.total_degree()) for f in zeta],
        "zeta_coefficients_in_prime_field": True,
        "zeta_basis_on_120_symmetric_products": [
            [int(c) for c in row] for row in kernel
        ],
        "F49_modulus_ascending": [1, 0, 1],
        "splitting_field_degree_over_F7": 12,
        "splitting_field_modulus_ascending": list(map(int, splitting.modulus())),
        "diagonal_negative_control_maps_to_all_zero_vector": True,
        "experiments": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).with_name("results") / "quotient-reconstruction.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
