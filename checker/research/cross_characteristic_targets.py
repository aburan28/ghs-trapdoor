"""Follow one guided F2197 route through all twelve theta reconstructions.

The arithmetic globals in quotient_reconstruction are rebound only in this
dedicated process. This is an additional exact toy case, not a general backend
API or a source-to-divisor evaluator.
"""

import argparse
import json
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from flint import fq_default_ctx, fq_default_poly_ctx

from . import quotient_reconstruction as qr
from .small_fields import frobenius_polynomial
from .weil_quotient_kernels import coordinates, frobenius_matrix, lagrangians, row_times

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def count_quartic(p, degree, quartic):
    k = fq_default_ctx(p, degree, var="v")
    R = fq_default_poly_ctx(k)
    T = R.gen()
    q = p**degree
    ff = {m: int(c.to_list()[0]) for m, c in quartic.items()}

    def root_count(f):
        if f.is_zero():
            return q
        if f.degree() == 0:
            return 0
        return (T.pow_mod(q, f) - T).gcd(f).degree()

    count = 0
    z = k.gen()
    for n in range(q):
        x = sum(k((n // p**j) % p) * z**j for j in range(degree))
        c = [k(0) for _ in range(5)]
        for (i, j, h), coefficient in ff.items():
            c[j] += k(coefficient) * x**i
        count += root_count(R(c))
    infinity = [k(0) for _ in range(5)]
    for (i, j, h), coefficient in ff.items():
        if h == 0:
            infinity[i] += k(coefficient)
    return count + root_count(R(infinity)) + int(ff.get((4, 0, 0), 0) == 0)


def run():
    path = Path(__file__).with_name("results") / "eight-isogeny-search.json"
    census = json.loads(path.read_text())
    field = next(row for row in census["fields"] if row["p"] == 13)
    route = next(
        row
        for row in field["routes"]
        if row["subgroup_order"] == 13**3 + 1 - row["trace"]
    )
    p = 13
    start = perf_counter()
    qr.K = fq_default_ctx(p, 18, var="s")
    qr.R = fq_default_poly_ctx(qr.K)
    qr.z = min(
        (root for root, multiplicity in qr.R(field["field_modulus_ascending"]).roots()),
        key=lambda value: tuple(qr.coefficients(value)),
    )

    def elt(n):
        return sum(qr.K((n // p**i) % p) * qr.z**i for i in range(3))

    A, B = elt(route["target_a"]), elt(route["target_b"])
    roots = sorted(
        (x for x, multiplicity in qr.R([-A * A, 12 * B, 6 * A, 0, 3]).roots()),
        key=lambda value: tuple(qr.coefficients(value)),
    )
    assert len(roots) == 4
    w = qr.EllipticWeilSet(A, B, roots[0], roots[1])
    weils = [w, w.conjugate(1), w.conjugate(2)]
    setup = perf_counter() - start
    # All source points are evaluated through this selected degree-three edge.
    aa, bb, u = elt(route["source_a"]), elt(route["source_b"]), elt(route["kernel_x"])
    v, ww = 6 * u * u + 2 * aa, 4 * (u**3 + aa * u + bb)
    assert aa - 5 * v == A and bb - 7 * (ww + u * v) == B

    def phi(P):
        if P is None:
            return None
        x, y = P
        return x + v / (x - u) + ww / (x - u) ** 2, y * (
            1 - v / (x - u) ** 2 - 2 * ww / (x - u) ** 3
        )

    generator = tuple(elt(n) for n in route["generator"])
    image_generator = phi(generator)
    current = None
    images = set()
    for _ in range(route["subgroup_order"]):
        image = phi(current)
        images.add(image)
        if image is not None:
            assert image[1] ** 2 == image[0] ** 3 + A * image[0] + B
        following = qr.add(aa, current, generator)
        assert phi(following) == qr.add(A, image, image_generator)
        current = following
    assert current is None and len(images) == route["subgroup_order"]
    action = frobenius_matrix(p, 3, [[2, 0], [0, 2]])
    targets = []
    for index, (rows, intersections) in enumerate(lagrangians(3, p % 3)):
        if any(intersections) or any(
            coordinates(row_times(row, action, 3), rows, 3) is None for row in rows
        ):
            continue
        before = perf_counter()
        seed, _theta, gradients, attempts = qr.choose_seed(weils, rows)
        curve = qr.gradient_quartic(gradients, base_prime=p)
        assert all(c.frobenius(1) == c for c in curve.values())
        assert qr.smooth(curve)
        for line in qr.rational_bitangents(gradients, base_prime=p).values():
            restricted = qr.restrict_line(curve, line)
            scalar, factors = restricted.factor()
            assert all(multiplicity % 2 == 0 for factor, multiplicity in factors)
            assert (4 - restricted.degree()) % 2 == 0 and scalar.is_square()
        counts = [count_quartic(p, n, curve) for n in (1, 2, 3)]
        polynomial = frobenius_polynomial(p, counts)
        expected = [p**3, 0, 0, -route["trace"], 0, 0, 1]
        opposite = [((-1) ** i) * c for i, c in enumerate(expected)]
        assert polynomial in (expected, opposite)
        targets.append(
            {
                "kernel": index,
                "kernel_basis_rows": rows,
                "seed": seed,
                "seed_attempts": attempts,
                "quartic_coefficients_ascending_monomials": [
                    [m, int(c.to_list()[0])] for m, c in sorted(curve.items())
                ],
                "counts_over_F13_F169_F2197": counts,
                "full_Frobenius_polynomial_ascending": polynomial,
                "matches_source": polynomial == expected,
                "Jacobian_order": sum(polynomial),
                "elapsed_seconds": perf_counter() - before,
            }
        )
    assert len(targets) == 12
    return {
        "experiment": 1,
        "prime": p,
        "selected_route": route,
        "source_field_modulus_ascending": field["field_modulus_ascending"],
        "working_field_modulus_ascending": list(map(int, qr.K.modulus())),
        "setup_seconds": setup,
        "elapsed_seconds": perf_counter() - start,
        "source_points_transferred_and_additions_checked": len(images),
        "targets": targets,
        "matching_targets": sum(t["matches_source"] for t in targets),
        "opposite_sign_targets": sum(not t["matches_source"] for t in targets),
        "input_sha256": sha256(path.read_bytes()).hexdigest(),
        "scope": {
            "all_twelve_kernels_of_selected_source_reconstructed": True,
            "all_28_bitangents_per_target_checked": True,
            "divisor_map_to_these_targets_evaluated": False,
            "new_cryptographic_family_established": False,
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
