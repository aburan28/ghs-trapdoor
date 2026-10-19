"""Sage 10.6: lift the two toy quotient maps to actual Jacobian divisor classes.

Run with sage -python. Inputs are produced by quotient_reconstruction's
--export-divisor-inputs command. Target inversion enumerates the order-373
Jacobian group, so this is a finite correctness experiment with O(r) cost,
not an efficient cryptographic point-transfer algorithm.
"""

import argparse
import itertools
import json
import sys
from hashlib import sha256
from pathlib import Path

from sage.all import GF, Curve, EllipticCurve, PolynomialRing, matrix, prod, vector
from sage.env import SAGE_VERSION

if not __debug__:
    raise RuntimeError("The research certificates require assertions; omit -O")


def normalize(values):
    lead = next(x for x in values if x != 0)
    return tuple(x / lead for x in values)


def verify_source(data, row, progress=True):
    Fp = GF(7)
    polynomial_ring = PolynomialRing(Fp, "t")
    K = GF(
        7**18,
        name="s",
        modulus=polynomial_ring(data["working_field_modulus_ascending"]),
    )
    s = K.gen()

    def elt(coefficients):
        return sum(K(a) * s**i for i, a in enumerate(coefficients))

    def log(message):
        if progress:
            print(
                f"source {row['source_a']},{row['source_b']}: {message}",
                file=sys.stderr,
                flush=True,
            )

    R = PolynomialRing(Fp, names=("X", "Y", "Z"))
    X, Y, Z = R.gens()
    original = sum(Fp(c) * X ** m[0] * Y ** m[1] * Z ** m[2] for m, c in row["quartic"])
    # A monic affine model with a transverse line at infinity keeps the
    # function-field ideal computations small over the large constant field.
    for a, b, c in itertools.product(range(7), repeat=3):
        moved = original(X + a * Y + b * Z, Y + c * Z, Z)
        binary = polynomial_ring(
            [moved.monomial_coefficient(Y**i * Z ** (4 - i)) for i in range(5)]
        )
        if (
            moved(0, 0, 1) != 0
            and binary.degree() == 4
            and binary.gcd(binary.derivative()).degree() == 0
        ):
            break
    else:
        raise ArithmeticError("No supported monic projective chart")
    C0 = Curve(moved)
    assert C0.is_smooth() and C0.genus() == 3
    C = C0.change_ring(K)
    F = C.function_field()
    coords = C.coordinate_functions()
    assert coords[0] == 1, "Unsupported affine chart chosen by Sage"
    raw_lines = {
        int(n): list(map(elt, coeffs)) for n, coeffs in row["bitangents"].items()
    }
    line_vectors = {
        n: (v[0], a * v[0] + v[1], b * v[0] + c * v[1] + v[2])
        for n, v in raw_lines.items()
    }
    lines = {
        n: sum(v * x for v, x in zip(coefficients, coords))
        for n, coefficients in line_vectors.items()
    }
    rational_points = C0.rational_points()
    O0 = next(
        P
        for P in rational_points
        if P[0] != 0
        and all(
            sum(v * K(x) for v, x in zip(line, P)) != 0
            for line in line_vectors.values()
        )
    )
    P0 = next(P for P in rational_points if P != O0 and P[0] != 0)
    Oplace, Pplace = C(list(O0)).place(), C(list(P0)).place()
    O = Oplace.divisor()
    # The projection center is outside C, hence Y/X has precisely this
    # hyperplane pole divisor in the selected chart.
    H = sum(-m * p.divisor() for p, m in coords[1].divisor().list() if m < 0)
    assert H.degree() == 4
    contacts = {}
    for n, line in lines.items():
        divisor = line.divisor() + H
        assert divisor.degree() == 4 and all(m % 2 == 0 for p, m in divisor.list())
        contacts[n] = sum((m // 2) * p.divisor() for p, m in divisor.list())
    # The three odd characteristics have xor zero, so their contact
    # divisors minus a canonical hyperplane represent the even delta.
    delta = contacts[9] + contacts[19] + contacts[26] - H
    initial = delta + O
    origin_effective = initial + initial.basis_function_space()[0].divisor()
    assert origin_effective.degree() == 3 and origin_effective.is_effective()
    # Replace delta by a linearly equivalent compact divisor. With this
    # Hess base divisor, defining_divisor() directly supplies D with
    # j(D)=[D-delta-O], avoiding a new Riemann-Roch solve per evaluation.
    delta = origin_effective - O
    J = C.jacobian(model="hess", base_div=origin_effective)
    G = J.group()
    generator = G.point(Pplace - Oplace)
    assert generator != G.zero() and 373 * generator == G.zero()
    phi = list((O + 2 * delta).basis_function_space())
    even = [1, 2, 4, 8, 16, 32]
    psi, norm_functions, differences = [], [], []
    for n in even:
        aa, bb = next((aa, aa ^ n) for aa in sorted(contacts) if aa ^ n in contacts)
        difference = contacts[aa] - contacts[bb]
        differences.append(difference)
        psi.append(list((difference + O + 2 * delta).basis_function_space()))
        norm_functions.append(lines[aa] / lines[bb])
    assert len(phi) == 3 and all(len(basis) == 3 for basis in psi)
    functions = phi + [fn for basis in psi for fn in basis] + norm_functions
    allowed = [O + 2 * delta] * 3 + [
        O + 2 * delta + d for d in differences for _ in range(3)
    ]
    exponents = [dict(d.list()) for d in allowed] + [
        dict((-2 * d).list()) for d in differences
    ]
    compiled = {
        fn: [
            (list(coefficient.numerator()), list(coefficient.denominator()))
            for coefficient in fn.list()
        ]
        for fn in functions
    }
    residue_generators = {}
    descent_matrices = {}

    def evaluations(divisor):
        places = divisor.list()
        assert sum(p.degree() * m for p, m in places) == 3
        if any(m != 1 for p, m in places):
            raise ArithmeticError(
                "Repeated support needs confluent evaluation; not implemented"
            )
        high = max((p for p, m in places), key=lambda p: p.degree())
        L, _, _ = high.residue_field()
        sigma = L(s) if L.has_coerce_map_from(K) else F(s).evaluate(high)
        residue_generators[high] = sigma
        converted = {}

        def lift(constant):
            if L == K:
                return constant
            if L.has_coerce_map_from(K):
                return L(constant)
            if constant not in converted:
                value = L(0)
                for coefficient in reversed(constant.polynomial().list()):
                    value = value * sigma + L(coefficient)
                converted[constant] = value
            return converted[constant]

        def evaluate_polynomial(coefficients, y):
            value = L(0)
            for coefficient in reversed(coefficients):
                value = value * y + lift(coefficient)
            return value

        values = []
        for place, multiplicity in places:
            y, zz = [fn.evaluate(place) for fn in coords[1:]]
            if place.degree() == 1 and high.degree() > 1:
                y, zz = lift(y), lift(zz)
            column = []
            for fn, powers in zip(functions, exponents):
                exponent = powers.get(place, 0)
                regular = (
                    fn if exponent == 0 else fn * place.local_uniformizer() ** exponent
                )
                if regular not in compiled:
                    compiled[regular] = [
                        (list(coefficient.numerator()), list(coefficient.denominator()))
                        for coefficient in regular.list()
                    ]
                value = L(0)
                try:
                    for numerator, denominator in reversed(compiled[regular]):
                        value = value * zz + evaluate_polynomial(
                            numerator, y
                        ) / evaluate_polynomial(denominator, y)
                except ZeroDivisionError:
                    value = regular.evaluate(place)
                    if place.degree() == 1 and high.degree() > 1:
                        value = lift(value)
                column.append(value)
            for j in range(place.degree()):
                values.append([v ** (K.cardinality() ** j) for v in column])
        return matrix(values), high, lift

    def descend(values, high):
        L = values[0].parent()
        if L == K:
            return tuple(values)

        def prime_vector(value):
            coefficients = list(value.polynomial())
            return vector(Fp, coefficients + [0] * (L.degree() - len(coefficients)))

        key = (L, residue_generators[high])
        if key not in descent_matrices:
            powers = [residue_generators[high] ** i for i in range(18)]
            descent_matrices[key] = matrix(Fp, [prime_vector(v) for v in powers])
        change = descent_matrices[key]
        out = []
        for value in values:
            coefficients = change.solve_left(prime_vector(value))
            out.append(
                sum(K(coefficient) * s**i for i, coefficient in enumerate(coefficients))
            )
        return tuple(out)

    initial_values, high, _lift = evaluations(origin_effective)
    denominator = initial_values.matrix_from_columns([0, 1, 2]).det() ** 2
    raw = []
    for i in range(6):
        det = initial_values.matrix_from_columns(
            [3 + 3 * i + j for j in range(3)]
        ).det()
        raw.append(
            det**2 * prod(initial_values[j, 21 + i] for j in range(3)) / denominator
        )
    raw = descend(raw, high)
    beta = [elt(row["theta_nulls"][n]) / value for n, value in zip(even, raw)]

    def jacobian_image(point):
        M, high, lift = evaluations(point.defining_divisor())
        values = [M.matrix_from_columns([0, 1, 2]).det() ** 2]
        for i in range(6):
            determinant = M.matrix_from_columns([3 + 3 * i + j for j in range(3)]).det()
            values.append(
                lift(beta[i]) * determinant**2 * prod(M[j, 21 + i] for j in range(3))
            )
        return normalize(descend(normalize(values), high))

    source_images = []
    for encoded in row["source_seven_coordinate_images_hex"]:
        coefficients = bytes.fromhex(encoded)
        assert len(coefficients) == 126 and all(c < 7 for c in coefficients)
        source_images.append(
            normalize([elt(coefficients[18 * i : 18 * i + 18]) for i in range(7)])
        )
    assert len(source_images) == 373 and len(set(source_images)) == 187
    target_images, target_points, inverse = [None] * 373, [None] * 373, {}
    current = G.zero()
    for n in range(187):
        image = jacobian_image(current)
        assert image in source_images and image not in inverse
        inverse[image] = [n] if n == 0 else [n, 373 - n]
        target_images[n], target_images[-n % 373] = image, image
        target_points[n], target_points[-n % 373] = current, -current
        if n % 20 == 0:
            log(f"{n}/187 Jacobian classes evaluated")
        current += generator
    assert len(inverse) == 187
    anchor = min(inverse[source_images[1]])
    lifted = []
    for n, image in enumerate(source_images):
        candidates = [
            j
            for j in inverse[image]
            if target_images[(j + anchor) % 373] == source_images[(n + 1) % 373]
        ]
        assert len(candidates) == 1
        lifted.append(candidates[0])
    # Scalar labels are used only to verify the constructed lift here.
    assert all(j == n * anchor % 373 for n, j in enumerate(lifted))
    # Every target class is rational over F7: use the base-field Hess group
    # for the exhaustive addition checks, and also check a scalar relation
    # in the extension-field group whose Weil functions were evaluated.
    Obase, Pbase = O0.place(), P0.place()
    Gbase = C0.jacobian(model="hess", base_div=3 * Obase.divisor()).group()
    base_generator = Gbase.point(Pbase - Obase)
    assert 373 * base_generator == Gbase.zero()
    base_points, current_base = [], Gbase.zero()
    for _ in range(373):
        base_points.append(current_base)
        current_base += base_generator
    assert current_base == Gbase.zero()
    for n in range(373):
        assert (
            base_points[lifted[n]] + base_points[anchor]
            == base_points[lifted[(n + 1) % 373]]
        )
    assert 37 * target_points[anchor] == target_points[lifted[37]]
    log("all 373 sign-resolved images and generator additions verified")
    return {
        "source_a": row["source_a"],
        "source_b": row["source_b"],
        "kernel_enumeration_index": 153,
        "target_chart_matrix": [[1, a, b], [0, 1, c], [0, 0, 1]],
        "target_base_point": list(map(int, O0)),
        "target_generator_point": list(map(int, P0)),
        "target_generator_class": "[P0-O0]",
        "source_generator_image_scalar_in_target_generator": int(anchor),
        "source_generator_image_divisor_over_F7": str(
            (anchor * base_generator).divisor()
        ),
        "source_scalar_37_image_divisor_over_F7": str(
            (37 * anchor * base_generator).divisor()
        ),
        "Jacobian_classes_evaluated_modulo_sign": 187,
        "source_points_lifted": 373,
        "generator_addition_checks": 373,
        "generator_addition_check_field": "F7",
        "scalar_37_relation_verified": True,
        "all_source_scalar_relations_verified": True,
        "inversion_method": "exhaustive enumeration of the tiny Jacobian group modulo sign",
        "inversion_cost_linear_in_subgroup_order": True,
        "scalable_divisor_inversion_implemented": False,
    }


def verify_conjugate_source(data, first_result, progress=True):
    """Reuse the certified lift through the explicit Galois/scaling isomorphism."""
    Fp = GF(7)
    PP = PolynomialRing(Fp, "t")
    K = GF(7**18, name="s", modulus=PP(data["working_field_modulus_ascending"]))
    s = K.gen()
    z = sum(K(c) * s**i for i, c in enumerate(data["census_z_embedding_coefficients"]))

    def encode(n):
        return K(n % 7) + K(n // 7 % 7) * z + K(n // 49) * z * z

    first, second = data["sources"]
    E1 = EllipticCurve(K, [encode(first["source_a"]), encode(first["source_b"])])
    E2 = EllipticCurve(K, [encode(second["source_a"]), encode(second["source_b"])])
    u = encode(78)
    assert u**4 * E1.a4() ** 7 == E2.a4() and u**6 * E1.a6() ** 7 == E2.a6()
    g1, g2 = (
        E1([encode(c) for c in first["source_generator"]]),
        E2([encode(c) for c in second["source_generator"]]),
    )
    assert 373 * g1 == E1(0) and 373 * g2 == E2(0)

    def canonical_sqrt(value):
        root = value.sqrt()

        def key(x):
            coefficients = list(map(int, x.polynomial()))
            return tuple(coefficients + [0] * (18 - len(coefficients)))

        return min((root, -root), key=key)

    bases = []
    for E, xx in ((E1, (19, 65)), (E2, (36, 47))):
        bases.append(
            [
                (
                    encode(x),
                    canonical_sqrt(encode(x) ** 3 + E.a4() * encode(x) + E.a6()),
                )
                for x in xx
            ]
        )
    T = matrix(GF(3), 6, 6)
    for i in range(3):
        j = (i - 1) % 3
        EE = EllipticCurve(K, [E2.a4() ** (7**j), E2.a6() ** (7**j)])
        target_basis = [EE([c ** (7**j) for c in P]) for P in bases[1]]
        for b, P in enumerate(bases[0]):
            xx, yy = [c ** (7**i) for c in P]
            image = EE([u ** (2 * 7**j) * xx, u ** (3 * 7**j) * yy])
            coefficients = next(
                (a, b)
                for a, b in itertools.product(range(3), repeat=2)
                if a * target_basis[0] + b * target_basis[1] == image
            )
            T[2 * i + b, 2 * j], T[2 * i + b, 2 * j + 1] = coefficients
    source_kernel = matrix(GF(3), first["kernel_basis_rows"])
    target_kernel = (source_kernel * T).echelon_form()
    expected_397 = matrix(
        GF(3), [[1, 0, 0, 1, 2, 1], [0, 1, 0, 2, 0, 1], [0, 0, 1, 1, 1, 1]]
    )
    assert target_kernel == expected_397

    R = PolynomialRing(Fp, names=("X", "Y", "Z"))
    X, Y, Z = R.gens()
    original = sum(
        Fp(c) * X ** m[0] * Y ** m[1] * Z ** m[2] for m, c in first["quartic"]
    )
    M = matrix(Fp, first_result["target_chart_matrix"])
    moved = original(
        *[sum(M[i, j] * v for j, v in enumerate((X, Y, Z))) for i in range(3)]
    )
    C = Curve(moved)
    O, P = (
        C(first_result["target_base_point"]).place(),
        C(first_result["target_generator_point"]).place(),
    )
    G = C.jacobian(model="hess", base_div=3 * O.divisor()).group()
    generator = G.point(P - O)
    anchor = first_result["source_generator_image_scalar_in_target_generator"]
    source_points, image_table = [], {}
    pp = E1(0)
    image_generator = anchor * generator
    image = G.zero()
    for n in range(373):
        expected = first["source_points_for_validation"][n]
        assert pp == (E1(0) if expected is None else E1([encode(c) for c in expected]))
        source_points.append(pp)
        image_table[pp] = image
        pp += g1
        image += image_generator
    assert pp == E1(0) and image == G.zero()
    qq = E2(0)
    images = []
    for n in range(373):
        expected = second["source_points_for_validation"][n]
        assert qq == (E2(0) if expected is None else E2([encode(c) for c in expected]))
        inverse = (
            E1(0) if qq == E2(0) else E1([(qq[0] / u**2) ** 49, (qq[1] / u**3) ** 49])
        )
        assert inverse == source_points[n]
        images.append(image_table[inverse])
        qq += g2
    for n in range(373):
        assert images[n] + images[1] == images[(n + 1) % 373]
    assert 37 * images[1] == images[37]
    if progress:
        print(
            "source 3,76: all 373 composed divisor images verified",
            file=sys.stderr,
            flush=True,
        )
    return {
        "source_a": second["source_a"],
        "source_b": second["source_b"],
        "kernel_enumeration_index": 397,
        "kernel_basis_rows": [list(map(int, row)) for row in target_kernel.rows()],
        "reuses_certified_first_source_Jacobian_target": True,
        "new_Jacobian_inversions_required": 0,
        "source_points_transported_and_lifted": 373,
        "generator_addition_checks": 373,
        "generator_addition_check_field": "F7",
        "scalar_37_relation_verified": True,
        "input_transport_to_first_curve": "(x,y) -> ((x/u^2)^49,(y/u^3)^49), u=78",
        "source_generator_image_divisor_over_F7": str(images[1].divisor()),
        "source_scalar_37_image_divisor_over_F7": str(images[37].divisor()),
        "scalable_divisor_inversion_implemented": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).with_name("results") / "quotient-divisor-inputs.json",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--source-index", type=int, choices=(0, 1))
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    if args.source_index is None:
        first_result = verify_source(data, data["sources"][0])
        results = [first_result, verify_conjugate_source(data, first_result)]
    else:
        results = [verify_source(data, data["sources"][args.source_index])]
    report = {
        "scope": {
            "toy_only": True,
            "Sage_version_used": SAGE_VERSION,
            "enumerative_divisor_inversion": True,
            "independent_Weil_restriction_examples": 1,
            "practical_attack_speedup": False,
        },
        "input_sha256": sha256(args.input.read_bytes()).hexdigest(),
        "sources": results,
    }
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
