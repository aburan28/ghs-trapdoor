"""Bounded Artin-Schreier presentation reduction and a DLP-descent control.

Equivalent equations are related by y->y+h(x). A subfield equation alone
does not put extension-field Jacobian points in the smaller rational group.
Only odd binary field degrees and polynomial single-pole covers are handled.
"""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from ghs_checker.field import BinaryField, is_irreducible

from .small_fields import factors, frobenius_polynomial

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def trim(polynomial):
    polynomial = list(polynomial)
    while polynomial and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def add_polynomials(a, b):
    return trim(
        [
            (a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0)
            for i in range(max(len(a), len(b)))
        ]
    )


def coboundary(k, h):
    result = [0] * max(1, 2 * len(h) - 1)
    for i, c in enumerate(h):
        result[2 * i] ^= k.square(c)
        result[i] ^= c
    return trim(result)


def reduce_polynomial(k, f):
    """Return g,h with f=g+h^2+h; nonconstant terms of g have odd degree."""
    if k.degree % 2 == 0:
        raise ValueError(
            "This experiment's canonical constant reduction requires odd field degree"
        )
    f = [k.element(c) for c in f]
    original = trim(f)
    g = f + [0] if not f else list(f)
    h = [0] * max(1, (len(g) + 1) // 2)
    for n in range(len(g) - 1, 0, -1):
        if n % 2 or not g[n]:
            continue
        root = k.frobenius(g[n], k.degree - 1)
        h[n // 2] ^= root
        g[n] = 0
        g[n // 2] ^= root
    constant = g[0]
    half_trace = 0
    for i in range((k.degree + 1) // 2):
        half_trace ^= k.frobenius(constant, 2 * i)
    h[0] ^= half_trace
    g[0] = k.absolute_trace(constant)
    g, h = trim(g), trim(h)
    if add_polynomials(g, coboundary(k, h)) != original:
        raise ArithmeticError("Artin-Schreier polynomial identity failed")
    return g, h


def evaluate(k, polynomial, x):
    value = 0
    for c in reversed(polynomial):
        value = k.mul(value, x) ^ c
    return value


def field(degree):
    modulus = next(
        m for m in range(2**degree + 1, 2 ** (degree + 1), 2) if is_irreducible(m)
    )
    return BinaryField(degree, modulus)


def embed_cubic(k, element):
    if k.degree % 3:
        raise ValueError("Need an extension of F8")
    root = next(z for z in range(k.size) if k.pow(z, 3) ^ z ^ 1 == 0)
    return (
        (element & 1)
        ^ k.mul((element >> 1) & 1, root)
        ^ k.mul((element >> 2) & 1, k.square(root))
    )


def curve_count(k, polynomial):
    # A nonconstant odd pole order gives exactly one smooth point at infinity.
    return 1 + 2 * sum(
        k.absolute_trace(evaluate(k, polynomial, x)) == 0 for x in range(k.size)
    )


def power_sums(polynomial, bound):
    degree = len(polynomial) - 1
    descending = list(reversed(polynomial))
    sums = [degree]
    for n in range(1, bound + 1):
        if n <= degree:
            sums.append(
                -sum(descending[j] * sums[n - j] for j in range(1, n))
                - n * descending[n]
            )
        else:
            sums.append(-sum(descending[j] * sums[n - j] for j in range(1, degree + 1)))
    return sums


def run():
    base = BinaryField(3, 0b1011)
    extension = field(9)
    selected = None
    # Fix a non-subfield coefficient family, and retain the first case with
    # a prime in the extension group absent from the base rational group.
    for a in range(2, 8):
        for b in range(1, 8):
            base_f = [0, b, 0, a, 0, 1]
            over6 = field(6)
            counts = [
                curve_count(base, base_f),
                curve_count(over6, [embed_cubic(over6, c) for c in base_f]),
            ]
            polynomial = frobenius_polynomial(8, counts)
            sums = power_sums(polynomial, 6)
            extension_poly = frobenius_polynomial(
                512, [513 - sums[3], 512**2 + 1 - sums[6]]
            )
            base_order, extension_order = sum(polynomial), sum(extension_poly)
            absent = [
                r for r, _ in factors(extension_order) if r > 31 and base_order % r
            ]
            if absent:
                selected = (
                    base_f,
                    counts,
                    polynomial,
                    extension_poly,
                    base_order,
                    extension_order,
                    absent,
                    sums,
                )
                break
        if selected:
            break
    if selected is None:
        raise ArithmeticError("No subgroup-size control in the bounded family")
    (
        base_f,
        counts,
        polynomial,
        extension_poly,
        base_order,
        extension_order,
        absent,
        sums,
    ) = selected
    normal = [embed_cubic(extension, c) for c in base_f]
    actual = curve_count(extension, normal)
    assert actual == 513 - sums[3]
    experiments = []
    for n in range(16):
        digest = sha256(f"ghs-artin-schreier-classes-v1:{n}".encode()).digest()
        h = [int.from_bytes(digest[2 * i : 2 * i + 2], "big") % 512 for i in range(7)]
        h[-1] = h[-1] or 1
        raw = add_polynomials(normal, coboundary(extension, h))
        reduced, correction = reduce_polynomial(extension, raw)
        expected, _ = reduce_polynomial(extension, normal)
        assert reduced == expected
        assert all(extension.frobenius(c, 3) == c for c in reduced)
        assert any(extension.frobenius(c, 3) != c for c in raw)
        # Check both y-values at every x using the polynomial identity. No
        # search for a square root is needed to verify the bijection.
        checked = 0
        square_image = {extension.square(y) ^ y: [] for y in range(512)}
        for y in range(512):
            square_image[extension.square(y) ^ y].append(y)
        for x in range(512):
            yy = evaluate(extension, raw, x)
            hh = evaluate(extension, correction, x)
            for y in square_image.get(yy, ()):
                Y = y ^ hh
                assert extension.square(Y) ^ Y == evaluate(extension, reduced, x)
                checked += 1
        assert checked + 1 == actual
        experiments.append(
            {
                "raw_polynomial_ascending": raw,
                "reduced_polynomial_ascending": reduced,
                "coordinate_correction_h_ascending": correction,
                "raw_polynomial_degree": len(raw) - 1,
                "reduced_odd_pole_order": len(reduced) - 1,
                "smooth_cover_genus": (len(reduced) - 2) // 2,
                "affine_point_transports_checked": checked,
                "polynomial_identity_verified": True,
                "reduced_coefficients_fixed_by_Frobenius_8": True,
            }
        )
    return {
        "scope": {
            "toy_only": True,
            "field_degrees_supported_for_reduction": "odd",
            "polynomial_single_pole_Artin_Schreier_covers_only": True,
            "nonpolynomial_covers_or_multiple_poles_implemented": False,
            "coordinate_isomorphism_is_not_a_Jacobian_DLP_descent": True,
            "norm_conorm_subgroup_map_evaluated": False,
            "novel_vulnerability_established": False,
        },
        "field_F512": extension.describe(),
        "field_F8": base.describe(),
        "base_curve_polynomial_ascending": base_f,
        "smooth_cover_genus": 2,
        "curve_counts_F8_F64_F512": counts + [actual],
        "F8_Frobenius_polynomial_ascending": polynomial,
        "F512_Frobenius_polynomial_ascending_from_cubed_roots": extension_poly,
        "F8_Jacobian_order": base_order,
        "F512_Jacobian_order": extension_order,
        "F512_Jacobian_order_factorization": factors(extension_order),
        "trace_zero_factor_order": extension_order // base_order,
        "supersingular_genus_2_toy_control": polynomial == [64, 0, 0, 0, 1],
        "extension_group_odd_primes_absent_from_base_Jacobian_order": absent,
        "DLP_descent_control": "Every homomorphism to the base rational Jacobian kills a subgroup of one of the displayed absent prime orders",
        "experiments": experiments,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
