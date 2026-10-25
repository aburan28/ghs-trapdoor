"""Sage 10.6 exact controls and bounded searches for experiments 1, 3, 6, 7, 8.

Run from checker/: sage -python research/eight_geometry_sage.py --experiment N --output FILE
"""

import argparse
import itertools
import json
import sys
from collections import Counter
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from sage.all import GF, ZZ, Curve, EllipticCurve, PolynomialRing, set_random_seed
from sage.env import SAGE_VERSION

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from research.small_fields import frobenius_polynomial, integer_poly_mul

ROOT = Path(__file__).resolve().parents[1]


if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")

def log(message):
    print(message, file=sys.stderr, flush=True)


def encode(x):
    k = x.parent()
    p = int(k.characteristic())
    return sum(int(c) * p**i for i, c in enumerate(x.polynomial()))


def decode(k, n):
    p, z = int(k.characteristic()), k.gen()
    return sum(k((n // p**i) % p) * z**i for i in range(k.degree()))


def hyper_count(k, coefficients):
    f = PolynomialRing(k, "x")(coefficients)
    total = 1 if f.degree() % 2 else (2 if f.leading_coefficient().is_square() else 0)
    for x in k:
        y2 = f(x)
        total += 1 if not y2 else 2 if y2.is_square() else 0
    return int(total)


def isogenies():
    from research.cubic_curve_census import counts_for_module

    reports = []
    for p in (7, 13):
        start = perf_counter()
        k = GF(p**3, "z")
        q = int(k.cardinality())
        nonsquare = k.multiplicative_generator()
        stats, records, routes = Counter(), sha256(), []
        direct_orbits, rescued_orbits, classes = set(), set(), set()
        for n in range(q):
            j = decode(k, n)
            if j in (k(0), k(1728)):
                continue
            base = EllipticCurve(j=j)
            for twist, e in enumerate((base, base.quadratic_twist(nonsquare))):
                stats["models"] += 1
                order = int(e.cardinality())
                t = q + 1 - order
                ordinary = t % p != 0
                roots = e.division_polynomial(3).roots(multiplicities=False)
                repeated = (t * t - 4 * q) % 3 == 0
                if repeated and len(roots) == 4:
                    ev = 2 * t % 3
                    module = ((ev, 0), (0, ev))
                else:
                    module = ((0, 1), (-q % 3, t % 3))
                _, mixed = counts_for_module(p, 3, module)
                selective = (
                    ordinary and (t * t - 4 * q) % 9 == 0 and (q + t + 1) % 9 == 0
                )
                orbit = tuple(sorted({encode(j ** (p**i)) for i in range(3)}))
                classes.add(t)
                stats["ordinary_models"] += ordinary
                stats["direct_mixed_models"] += bool(mixed)
                stats["selective_models"] += selective
                if mixed:
                    direct_orbits.add((t, orbit))
                row = [n, twist, order, len(roots), module, mixed, selective]
                records.update((json.dumps(row, separators=(",", ":")) + "\n").encode())
                # All rational first steps are examined for the selected negative models.
                if not selective or mixed:
                    continue
                for u in sorted(roots, key=encode):
                    phi = e.isogeny(e.division_polynomial(3).parent().gen() - u)
                    target = phi.codomain()
                    target_roots = target.division_polynomial(3).roots(
                        multiplicities=False
                    )
                    stats["selected_negative_edges"] += 1
                    if len(target_roots) != 4:
                        continue
                    stats["rescued_models"] += 1
                    rescued_orbits.add((t, orbit))
                    primes = [int(r) for r, _ in ZZ(order).factor() if r > 3]
                    if not primes:
                        continue
                    r = max(primes)
                    generator = None
                    for g in e.gens():
                        candidate = (order // r) * g
                        if candidate:
                            generator = candidate
                            break
                    if generator is None:
                        for g in e.points():
                            candidate = (order // r) * g
                            if candidate:
                                generator = candidate
                                break
                    assert generator and r * generator == e(0)
                    image = phi(generator)
                    assert image and r * image == target(0)
                    samples = sorted(
                        {0, 1, 2, r - 1}
                        | {
                            int(ZZ(sha256(f"{p}:{n}:{i}".encode()).hexdigest(), 16) % r)
                            for i in range(12)
                        }
                    )
                    for scalar in samples:
                        assert phi(scalar * generator) == scalar * image
                        assert (
                            phi(scalar * generator + generator)
                            == scalar * image + image
                        )
                    routes.append(
                        {
                            "j": n,
                            "twist": twist,
                            "trace": t,
                            "source_a": encode(e.a4()),
                            "source_b": encode(e.a6()),
                            "kernel_x": encode(u),
                            "target_j": encode(target.j_invariant()),
                            "target_a": encode(target.a4()),
                            "target_b": encode(target.a6()),
                            "degree": 3,
                            "subgroup_order": r,
                            "generator": [encode(x) for x in generator[:2]],
                            "image": [encode(x) for x in image[:2]],
                            "scalar_controls": samples,
                        }
                    )
        reports.append(
            {
                "p": p,
                "field_modulus_ascending": list(map(int, k.modulus())),
                "quadratic_twist_parameter": encode(nonsquare),
                "counts": dict(stats),
                "ordinary_and_supersingular_trace_classes_combined": sorted(classes),
                "direct_Frobenius_j_orbits_by_trace": len(direct_orbits),
                "rescued_Frobenius_j_orbits_by_trace": len(rescued_orbits),
                "ordered_model_data_sha256": records.hexdigest(),
                "routes": routes,
                "elapsed_seconds": perf_counter() - start,
            }
        )
        log(
            f"experiment 1: F{q}: {stats['models']} models, {len(routes)} certified routes"
        )
    return {
        "fields": reports,
        "scope": {
            "all_declared_generic_j_and_twists_exhausted": True,
            "selected_negative_rational_first_steps_exhausted": True,
            "new_theta_targets_reconstructed": False,
            "least_total_cost_target_certified": False,
            "comparison": "selection yield; no blind-search speedup inferred from this census",
        },
    }


def quartic_curve(kernel, field):
    report = json.loads(
        (ROOT / "research/results/quotient-reconstruction.json").read_text()
    )
    target = next(
        t
        for t in report["sources"][0]["targets"]
        if t["kernel_enumeration_index"] == kernel
    )
    ring = PolynomialRing(GF(7), names=("X", "Y", "Z"))
    X, Y, Z = ring.gens()
    original = sum(
        GF(7)(c) * X ** m[0] * Y ** m[1] * Z ** m[2]
        for m, c in target["quartic_coefficients_ascending_monomials"]
    )
    univ = PolynomialRing(GF(7), "t")
    for a, b, c in itertools.product(range(7), repeat=3):
        moved = original(X + a * Y + b * Z, Y + c * Z, Z)
        binary = univ(
            [moved.monomial_coefficient(Y**i * Z ** (4 - i)) for i in range(5)]
        )
        if (
            moved(0, 0, 1)
            and binary.degree() == 4
            and binary.gcd(binary.derivative()).degree() == 0
        ):
            break
    return Curve(moved).change_ring(field), target, [a, b, c]


def lift_signs():
    k = GF(49, "w", modulus=PolynomialRing(GF(7), "t")([1, 0, 1]))
    reports = []
    for kernel in (153, 810):
        C, target, chart = quartic_curve(kernel, k)
        points = C.rational_points()
        rational = [P for P in points if all(x**7 == x for x in P)]
        O = rational[0].place()
        G = C.jacobian(model="hess", base_div=3 * O.divisor()).group()
        polynomial = target["frobenius_polynomial_ascending"]
        plus, minus = (
            sum(polynomial),
            sum((-1) ** i * c for i, c in enumerate(polynomial)),
        )
        checks = []
        for P in points:
            if P in rational:
                D = G.point(P.place() - O)
                assert plus * D == G.zero()
                sigma = G.point(C([x**7 for x in P]).place() - O)
                assert sigma == D
                checks.append(
                    {"P": [encode(x) for x in P], "sign": 1, "annihilator": plus}
                )
            else:
                Q = C([x**7 for x in P])
                D = G.point(P.place() - Q.place())
                sigma = G.point(Q.place() - P.place())
                assert D != G.zero() and sigma == -D
                assert minus * D == G.zero()
                order = minus
                for prime, _ in ZZ(minus).factor():
                    while order % prime == 0 and (order // prime) * D == G.zero():
                        order //= int(prime)
                assert order > 1 and order * D == G.zero()
                checks.append(
                    {
                        "P": [encode(x) for x in P],
                        "Frob_P": [encode(x) for x in Q],
                        "sign": -1,
                        "order": order,
                        "divisor": str(D.divisor()),
                        "unordered_pair_of_classes_Frobenius_stable": True,
                    }
                )
        reports.append(
            {
                "kernel": kernel,
                "chart": chart,
                "invariant_order": plus,
                "anti_invariant_order": minus,
                "curve_points": len(points),
                "rational_points": len(rational),
                "checks": checks,
            }
        )
        log(f"experiment 3: kernel {kernel}: all {len(points)} point controls verified")
    return {
        "targets": reports,
        "scope": {
            "all_F49_curve_point_controls_exhausted": True,
            "actual_Jacobian_lift_signs_checked": True,
            "intrinsic_Kummer_descent_certified_by_unordered_class_pair": True,
            "explicit_F7_Kummer_coordinate_embedding_computed": False,
            "global_quotient_map_lift_sign_certificate_constructed": False,
            "source_to_opposite_sign_target_map_constructed": False,
        },
    }


def prym():
    cases = []
    p = 3
    for a, b in itertools.product(range(p), range(1, p)):
        start = perf_counter()
        f = [b, a] + [0] * 8 + [1]
        F = PolynomialRing(GF(p), "u")(f)
        assert F.gcd(F.derivative()).degree() == 0
        cover = [0] * 21
        cover[::2] = f
        quotient = [0] + f
        counts_c, counts_b, counts_d = [], [], []
        for n in range(1, 10):
            k = GF(p**n, "z") if n > 1 else GF(p)
            nc = hyper_count(k, cover)
            nb, nd = hyper_count(k, quotient), hyper_count(k, f)
            assert nc == nb + nd - int(k.cardinality()) - 1
            counts_c.append(nc)
            counts_b.append(nb)
            counts_d.append(nd)
        pc = frobenius_polynomial(p, counts_c)
        pb = frobenius_polynomial(p, counts_b[:5])
        pd = frobenius_polynomial(p, counts_d[:4])
        assert pc == integer_poly_mul(pb, pd)
        assert sum(pc) == sum(pb) * sum(pd)
        # Finite affine maps and free involution control over F9.
        k = GF(9, "z")
        checked = 0
        ff = PolynomialRing(k, "x")(f)
        for x in k:
            y2 = ff(x * x)
            for y in y2.sqrt(all=True) if y2.is_square() else []:
                assert y * y == ff(x * x) and (x * y) ** 2 == x * x * ff(x * x)
                assert not (x == -x and y == -y)
                checked += 1
        cases.append(
            {
                "p": p,
                "f_ascending": f,
                "C_rhs_ascending": cover,
                "B_rhs_ascending": quotient,
                "genera_C_B_D": [9, 5, 4],
                "counts_C_over_p_through_p9": counts_c,
                "counts_B_over_p_through_p9": counts_b,
                "counts_D_over_p_through_p9": counts_d,
                "full_Frobenius_polynomials_ascending": {"C": pc, "B": pb, "D": pd},
                "orders": {"C": sum(pc), "B": sum(pb), "Prym": sum(pd)},
                "Prym_order_factorization": [
                    [int(r), int(e)] for r, e in ZZ(sum(pd)).factor()
                ],
                "affine_point_maps_checked_F9": checked,
                "free_involution": "(x,y)->(-x,-y); f(0)!=0 excludes affine fixed points; infinity points exchanged",
                "Prym_principal_polarization": "etale double cover; restriction of theta_C is twice the canonical principal polarization",
                "Prym_Jacobian_identification": "pullback Jac(D)->Prym(C/B); C->D is ramified, hence pullback injective",
                "odd_subgroup_preservation": "pushforward_D o pullback_D = [2]",
                "elapsed_seconds": perf_counter() - start,
            }
        )
        log(f"experiment 7: a={a}, b={b}: all nine extension counts verified")
    return {
        "cases": cases,
        "bounds": {"prime": 3, "a": [0, 1, 2], "b": [1, 2]},
        "scope": {
            "all_declared_covers_exhausted": True,
            "full_zeta_factorization_independently_counted": True,
            "general_restriction_to_Prym_matching_search_completed": False,
            "pointwise_Prym_Jacobian_transfer_computed": False,
            "new_relation_algorithm_implemented": False,
        },
    }


def dimension_four():
    # The existing explicit cover is the constructive Jacobian-locus control.
    data = json.loads((ROOT / "research/results/toy-cover.json").read_text())
    F7 = GF(7)
    k = GF(49, "i", modulus=PolynomialRing(F7, "t")([1, 0, 1]))
    i = k.gen()
    ring = PolynomialRing(k, names=("x", "y"))
    x, y = ring.gens()
    h = sum(k(c) * x**j for j, c in enumerate(data["cover"]["rhs"]))
    f = (x - (1 + i)) * x * (x - 1) * (x - 2) * (x - 3)
    H, C = (
        Curve(y * y - h).projective_closure(2),
        Curve(y * y - f).projective_closure(2),
    )
    assert H.genus() == 4 and C.genus() == 2
    # Use finite split fibers to avoid concealing exceptional-chart work.
    source_points = []
    for P in C.rational_points():
        if not P[2] or not P[1]:
            continue
        xx, yy = P[0] / P[2], P[1] / P[2]
        univ = PolynomialRing(k, "t")
        t = univ.gen()
        roots = (t * t + 2 * (xx - 1) * t - 1).roots(multiplicities=False)
        if len(roots) != 2 or any(not tt or not 1 - i * tt for tt in roots):
            continue
        fibers = [H([tt, 8 * tt**3 * yy / (1 - i * tt), 1]) for tt in roots]
        assert all((-(Q[0] ** 2) + 2 * Q[0] + 1) / (2 * Q[0]) == xx for Q in fibers)
        source_points.append((P, fibers))
    assert len(source_points) >= 2
    O, Ofiber = source_points[0]
    Hbase = next(Q for Q in H.rational_points() if Q[2]).place()
    GC = C.jacobian(model="hess", base_div=2 * O.place().divisor()).group()
    GH = H.jacobian(model="hess", base_div=4 * Hbase.divisor()).group()
    n = int(data["jacobian_group_order"])
    witness = None
    for P, fiber in source_points[1:]:
        D = (n // 17) * GC.point(P.place() - O.place())
        pull = sum(Q.place().divisor() for Q in fiber) - sum(
            Q.place().divisor() for Q in Ofiber
        )
        conjugate = sum(H([v**7 for v in Q]).place().divisor() for Q in fiber) - sum(
            H([v**7 for v in Q]).place().divisor() for Q in Ofiber
        )
        image = (n // 17) * GH.point(pull + conjugate)
        if D != GC.zero() and image != GH.zero():
            assert 17 * D == GC.zero() and 17 * image == GH.zero()
            witness = {
                "P": [encode(v) for v in P],
                "O": [encode(v) for v in O],
                "cofactor": n // 17,
                "order": 17,
                "source_divisor": str(D.divisor()),
                "image_divisor": str(image.divisor()),
                "map": "D -> pi^*D + Frob7(pi^*D)",
                "finite_split_fiber_points": [[encode(v) for v in Q] for Q in fiber],
            }
            break
    assert witness is not None
    return {
        "constructive_genus4_control": data["frobenius_polynomials_ascending"],
        "order17_subgroup_witness": witness,
        "finite_split_fiber_source_points": len(source_points),
        "scope": {
            "actual_genus4_Jacobian_subgroup_image_verified": True,
            "general_Schottky_filter_implemented": False,
            "general_dimension4_quotient_reconstruction_implemented": False,
            "exceptional_and_nonsplit_map_charts_evaluated": False,
        },
    }


def wild_control():
    # A rank-two Frobenius-stable coefficient space over F4, embedded in F16.
    k = GF(16, "z")
    a = next(v for v in k if v not in (0, 1) and v**4 == v)
    assert a * a == a + 1
    controls = []
    for n in (1, 2, 3):
        L = GF(16**n, "s") if n > 1 else k
        embed = (
            k.hom([next(v for v in L if k.modulus()(v) == 0)], L)
            if n > 1
            else k.hom([k.gen()], k)
        )
        aa = embed(a)

        def trace(v, degree=4 * n):
            return sum(v ** (2**j) for j in range(degree))

        counts = [1, 1, 1, 1]
        for x in L:
            bits = [trace(c * x**3) == 0 for c in (aa, aa * aa, aa + aa * aa)]
            counts[0] += 4 * int(bits[0] and bits[1])
            for j in range(3):
                counts[j + 1] += 2 * int(bits[j])
        assert counts[0] == sum(counts[1:]) - 2 * (16**n + 1)
        controls.append(counts)
    ph = frobenius_polynomial(16, [c[0] for c in controls])
    factors = [[16, controls[0][j] - 16 - 1, 1] for j in (1, 2, 3)]
    expected = integer_poly_mul(integer_poly_mul(factors[0], factors[1]), factors[2])
    assert ph == expected
    return {
        "F16_modulus_ascending": list(map(int, k.modulus())),
        "a": encode(a),
        "equations": ["y1^2+y1=a*x^3", "y2^2+y2=a^2*x^3"],
        "Frobenius2_action": "a->a^2, y1->y2, y2->y1; a+a^2=1",
        "genus": 3,
        "quotient_genera": [1, 1, 1],
        "counts_H_C1_C2_C3_over_F16_F256_F4096": controls,
        "full_Frobenius_polynomial_H_ascending": ph,
        "quotient_Frobenius_polynomials_ascending": factors,
        "norm_conorm_identity": "pi_* pi^*=[2] for each degree-two quotient",
        "scope": {
            "explicit_wild_cover_equations": True,
            "all_three_full_quotient_polynomials_verified": True,
            "full_genus3_polynomial_independently_counted": True,
            "Jacobian_subgroup_norm_conorm_evaluated": False,
            "ordinary_binary_elliptic_curve_example": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--experiment", type=int, choices=(1, 3, 6, 7, 8), required=True
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    set_random_seed(20261020)
    start = perf_counter()
    functions = {
        1: isogenies,
        3: lift_signs,
        6: dimension_four,
        7: prym,
        8: wild_control,
    }
    result = functions[args.experiment]()
    result.update(
        {
            "experiment": args.experiment,
            "sage_version": SAGE_VERSION,
            "random_seed": 20261020,
            "elapsed_seconds": perf_counter() - start,
            "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "input_sha256": {
                str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                for p in (
                    ROOT / "research/results/quotient-reconstruction.json",
                    ROOT / "research/results/toy-cover.json",
                )
            },
            "practical_speedup_established": False,
        }
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
