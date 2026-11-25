"""Independent Sage field/elliptic oracles for the native Artin--Schreier sweep."""

import argparse
import itertools
import json
import sys
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from sage.all import GF, EllipticCurve, PolynomialRing, matrix
from sage.env import SAGE_VERSION

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ghs_checker.field import BinaryField
from research.wild_artin_schreier import class_span

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def dimension_four_controls():
    root = Path(__file__).with_name("results")
    saved = json.loads((root / "eight-dimension-four-kernels.json").read_text())
    F2 = GF(2)
    checked = 0
    for module in saved["modules"]:
        blocks = [2, 2, 2, 2] if module["module"].startswith("elliptic") else [4, 4]
        H = matrix(F2, 8, 8)
        offset = 0
        for d in blocks:
            for i, j in itertools.product(range(d), repeat=2):
                H[offset + i, offset + j] = int(i != j)
            offset += d
        F = matrix(
            F2,
            [[(v >> j) & 1 for j in range(8)] for v in module["Frobenius_rows_binary"]],
        )
        assert F * H * F.transpose() == H
        for kernel in module["kernels"]:
            M = matrix(
                F2,
                [[(v >> j) & 1 for j in range(8)] for v in kernel["basis_rows_binary"]],
            )
            assert M.rank() == 4 and M * H * M.transpose() == 0
            assert M.stack(M * F).rank() == 4
            offset = 0
            dimensions = []
            for d in blocks:
                columns = list(range(offset)) + list(range(offset + d, 8))
                dimensions.append(4 - M.matrix_from_columns(columns).rank())
                offset += d
            assert dimensions == kernel["factor_intersection_dimensions"]
            checked += 1
    # Realize each GL2(F2) conjugacy type on an actual ordinary E/F625.
    k = GF(625, "w")
    realized = {}
    for j in k:
        if j**5 == j or j in (k(0), k(1728)):
            continue
        E = EllipticCurve(j=j)
        if not E.is_ordinary():
            continue
        roots = E.division_polynomial(2).roots(multiplicities=False)
        kind = {0: "irreducible_order3", 1: "unipotent", 3: "identity"}[len(roots)]
        if kind not in realized:
            realized[kind] = {
                "j_coefficients": list(map(int, j.polynomial())),
                "a_coefficients": list(map(int, E.a4().polynomial())),
                "b_coefficients": list(map(int, E.a6().polynomial())),
                "order": int(E.cardinality()),
                "rational_2division_roots": len(roots),
            }
        if len(realized) == 3:
            break
    assert len(realized) == 3
    genus2_realized = []
    L = GF(25, "v")
    polynomials = PolynomialRing(L, "x")
    for module in saved["modules"]:
        partition = module.get("branch_cycle_partition")
        if partition is None:
            continue
        factors = []
        for degree in partition:
            for attempt in range(100):
                factor = polynomials.irreducible_element(degree)
                if factor not in factors:
                    factors.append(factor)
                    break
            else:
                raise ArithmeticError("Could not generate distinct irreducible factors")
        f = polynomials.one()
        for factor in factors:
            f *= factor
        assert f.degree() == 6 and f.gcd(f.derivative()).degree() == 0
        assert (
            sorted(factor.degree() for factor, multiplicity in f.factor()) == partition
        )
        genus2_realized.append(
            {
                "cycle_partition": partition,
                "F25_modulus_ascending": list(map(int, L.modulus())),
                "rhs_coefficients_in_polynomial_basis": [
                    list(map(int, c.polynomial())) for c in f
                ],
                "equation": "y^2=f(x)",
                "genus": 2,
                "branch_Frobenius_cycle_type_verified_by_factorization": True,
            }
        )
    assert len(genus2_realized) == 11
    return {
        "stored_kernels_independently_verified": checked,
        "F625_modulus_ascending": list(map(int, k.modulus())),
        "elliptic_module_realizations": realized,
        "genus2_all_Sp4_F2_cycle_type_realizations": genus2_realized,
        "genus2_module_realization": "the five finite roots and infinity of the saved weak F49 model are all rational, hence E[2] action I4",
        "Schottky_locus_not_inferred_from_kernel_checks": True,
    }


def cross_characteristic_controls():
    path = (
        Path(__file__).with_name("results") / "eight-cross-characteristic-targets.json"
    )
    data = json.loads(path.read_text())
    route = data["selected_route"]
    p = 13
    k = GF(
        p**3,
        "z",
        modulus=PolynomialRing(GF(p), "t")(data["source_field_modulus_ascending"]),
    )

    def decode(n):
        return sum(k((n // p**i) % p) * k.gen() ** i for i in range(3))

    E = EllipticCurve(k, [decode(route["source_a"]), decode(route["source_b"])])
    R = PolynomialRing(k, "x")
    phi = E.isogeny(R.gen() - decode(route["kernel_x"]))
    assert E.cardinality() == 2137
    assert phi.codomain().a4() == decode(
        route["target_a"]
    ) and phi.codomain().a6() == decode(route["target_b"])
    G = E([decode(n) for n in route["generator"]])
    image_G = phi(G)
    P = E(0)
    for _ in range(2137):
        assert phi(P + G) == phi(P) + image_G
        P += G
    assert P == E(0)
    count_checks = []
    for target in data["targets"]:
        counts = []
        for degree in (1, 2, 3):
            L = GF(p**degree, "v") if degree > 1 else GF(p)
            S = PolynomialRing(L, "T")
            count = 0
            ff = target["quartic_coefficients_ascending_monomials"]

            def roots(poly, field_size=p**degree):
                return (
                    field_size if poly == 0 else len(poly.roots(multiplicities=False))
                )

            for x in L:
                coefficients = [L(0)] * 5
                for (i, j, h), c in ff:
                    coefficients[j] += L(c) * x**i
                count += roots(S(coefficients))
            infinity = [L(0)] * 5
            for (i, j, h), c in ff:
                if h == 0:
                    infinity[i] += L(c)
            count += roots(S(infinity))
            count += int(not any(m == [4, 0, 0] and c != 0 for m, c in ff))
            counts.append(int(count))
        assert counts == target["counts_over_F13_F169_F2197"]
        count_checks.append({"kernel": target["kernel"], "independent_counts": counts})
    return {
        "all_2137_isogeny_additions_independently_verified": True,
        "all_twelve_targets_three_counts_independently_verified": count_checks,
        "input_sha256": sha256(path.read_bytes()).hexdigest(),
    }


def run():
    rows = []
    for n, modulus in ((4, 19), (6, 67)):
        native = BinaryField(n, modulus)
        F2 = GF(2)
        R = PolynomialRing(F2, "t")
        k = GF(2**n, "z", modulus=R([(modulus >> i) & 1 for i in range(n + 1)]))
        z = k.gen()
        values = [sum(k((v >> i) & 1) * z**i for i in range(n)) for v in range(2**n)]
        for a, b in itertools.product(range(2**n), repeat=2):
            assert values[native.mul(a, b)] == values[a] * values[b]
        trace = [native.absolute_trace(v) for v in range(2**n)]
        constant = next(v for v in range(2**n) if trace[v])
        tested, digest = 0, sha256()
        # Exhaust F16; cover every reported span type over F64.
        if n == 4:
            inputs = list(itertools.product(range(16), range(16), (0, 1)))
        else:
            saved = json.loads(
                (
                    Path(__file__).with_name("results") / "eight-wild-classes.json"
                ).read_text()
            )
            field = next(f for f in saved["fields"] if f["degree"] == n)
            inputs = sorted(
                {
                    (r["a"], r["b"], r["constant_trace_bit"])
                    for family in field["families"]
                    for r in family["representatives"]
                }
            )
        for a, b, bit in inputs:
            aa, bb, cc = values[a], values[b], values[constant] if bit else k(0)
            basis = class_span(native, a, b, bit)
            orbit = []
            av, bv = aa, bb
            for _ in range(n):
                orbit.append(
                    [av.polynomial()[i] for i in range(n)]
                    + [bv.polynomial()[i] for i in range(n)]
                    + [F2(bit)]
                )
                av, bv = av * av, bv * bv
            assert matrix(F2, orbit).rank() == len(basis)
            if a:
                E = EllipticCurve(k, [0, 0, aa, aa * bb, aa * aa * cc])
                count = 1 + 2 * sum(
                    trace[native.mul(a, native.pow(x, 3)) ^ native.mul(b, x)] == bit
                    for x in range(2**n)
                )
                assert E.cardinality() == count
                tested += 1
                digest.update(f"cubic:{a}:{b}:{bit}:{count}\n".encode())
            if a and b:
                E = EllipticCurve(k, [1, cc, 0, aa * bb, 0])
                assert E.is_ordinary()
                count = 2 + 2 * sum(
                    trace[native.mul(a, x) ^ native.mul(b, native.inv(x))] == bit
                    for x in range(1, 2**n)
                )
                assert E.cardinality() == count
                tested += 1
                digest.update(f"two_poles:{a}:{b}:{bit}:{count}\n".encode())
        rows.append(
            {
                "degree": n,
                "modulus_binary": modulus,
                "multiplication_pairs_verified": (2**n) ** 2,
                "class_inputs_checked": len(inputs),
                "genus1_point_counts_verified": tested,
                "ordered_count_sha256": digest.hexdigest(),
                "selection": "all classes"
                if n == 4
                else "all reported span-type representatives",
            }
        )
    return {
        "checks": rows,
        "dimension_four_controls": dimension_four_controls(),
        "cross_characteristic_controls": cross_characteristic_controls(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    start = perf_counter()
    result = run()
    result.update(
        {
            "sage_version": SAGE_VERSION,
            "elapsed_seconds": perf_counter() - start,
            "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "input_sha256": sha256(
                (
                    Path(__file__).with_name("results") / "eight-wild-classes.json"
                ).read_bytes()
            ).hexdigest(),
        }
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
