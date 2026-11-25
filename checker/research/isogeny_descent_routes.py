"""Exact small 3-isogeny graph and routing to an existing descent example.

This tests a construction-specific false negative, not a new attack. All
arithmetic is over F343; no cryptographic parameter search is performed.
"""

import argparse
import json
from collections import Counter, deque
from hashlib import sha256
from pathlib import Path

from .cubic_curve_census import counts_for_module, curve_data
from .small_fields import SmallField, ec_add, ec_mul, points


def quotient(k, a, b, u):
    """Normalized Velu codomain for kernel {O,+/-Q}, x(Q)=u, ord(Q)=3."""
    psi3 = k.sub(
        k.add(
            k.add(k.mul(3, k.pow(u, 4)), k.mul(6 % k.p, k.mul(a, k.pow(u, 2)))),
            k.mul(12 % k.p, k.mul(b, u)),
        ),
        k.mul(a, a),
    )
    if psi3:
        raise ValueError("u is not a root of the 3-division polynomial")
    v = k.add(k.mul(6 % k.p, k.pow(u, 2)), k.mul(2, a))
    w = k.mul(4, k.add(k.add(k.pow(u, 3), k.mul(a, u)), b))
    return k.sub(a, k.mul(5 % k.p, v)), k.sub(b, k.mul(7 % k.p, k.add(w, k.mul(u, v))))


def evaluate_isogeny(k, a, b, u, point, scale=1):
    """Velu map followed by (x,y)->(scale^2*x,scale^3*y).

    The kernel is rational as a subgroup; its nonidentity points need not be
    rational. These rational formulas avoid adjoining their y-coordinates.
    """
    if point is None or point[0] == u:
        return None
    x, y = point
    d_inv = k.inv(k.sub(x, u))
    v = k.add(k.mul(6 % k.p, k.pow(u, 2)), k.mul(2, a))
    w = k.mul(4, k.add(k.add(k.pow(u, 3), k.mul(a, u)), b))
    xx = k.add(x, k.add(k.mul(v, d_inv), k.mul(w, k.pow(d_inv, 2))))
    yy = k.mul(
        y,
        k.sub(k.sub(1, k.mul(v, k.pow(d_inv, 2))), k.mul(2, k.mul(w, k.pow(d_inv, 3)))),
    )
    return k.mul(k.pow(scale, 2), xx), k.mul(k.pow(scale, 3), yy)


def scaling(k, source, target):
    """Choose a stated K-isomorphism, rather than deduplicating by j alone."""
    for u in range(1, k.size):
        if (k.mul(k.pow(u, 4), source[0]), k.mul(k.pow(u, 6), source[1])) == target:
            return u
    raise ArithmeticError("Same-j models did not admit the required K-isomorphism")


def on_curve(k, a, b, point):
    return point is None or k.pow(point[1], 2) == k.add(
        k.add(k.pow(point[0], 3), k.mul(a, point[0])), b
    )


def run():
    k = SmallField(7, 3, [1, 0, 1, 1])
    k.tables()
    first = curve_data(k, 248, 42)
    nodes, edges, pending = {first["j"]: first}, [], deque([first["j"]])
    while pending:
        j = pending.popleft()
        source = nodes[j]
        a, b = source["a"], source["b"]
        for u in source["three_division_roots"]:
            aa, bb = quotient(k, a, b, u)
            target = curve_data(k, aa, bb)
            jj = target["j"]
            if jj not in nodes:
                nodes[jj] = target
                pending.append(jj)
            if len(nodes) > 32:
                raise ArithmeticError("Toy graph bound exceeded")
            target = nodes[jj]
            scale = scaling(k, (aa, bb), (target["a"], target["b"]))
            edges.append(
                {
                    "source_j": j,
                    "target_j": jj,
                    "kernel_x": u,
                    "normalized_codomain_a": aa,
                    "normalized_codomain_b": bb,
                    "codomain_scaling": scale,
                }
            )

    # Check every rational point and a generating addition on every edge.
    for edge in edges:
        source, target = nodes[edge["source_j"]], nodes[edge["target_j"]]
        a, b = source["a"], source["b"]
        group = [None] + list(points(k, [b, a, 0, 1]))
        generator = next(p for p in group if p is not None)
        if source["group_order"] != 373 or ec_mul(k, 0, a, generator, 373) is not None:
            raise ArithmeticError("Expected prime-order source group")

        def phi(point, a=a, b=b, edge=edge):
            return evaluate_isogeny(
                k, a, b, edge["kernel_x"], point, edge["codomain_scaling"]
            )

        image_generator = phi(generator)
        images = [phi(p) for p in group]
        if len(set(images)) != 373 or not all(
            on_curve(k, target["a"], target["b"], p) for p in images
        ):
            raise ArithmeticError("Isogeny failed curve membership or injectivity")
        for point in group:
            if phi(ec_add(k, 0, a, point, generator)) != ec_add(
                k, 0, target["a"], phi(point), image_generator
            ):
                raise ArithmeticError("Isogeny failed a generator-addition identity")
        digest = sha256(json.dumps(images, separators=(",", ":")).encode()).hexdigest()
        edge.update(
            {
                "source_generator": generator,
                "image_generator": image_generator,
                "rational_point_images_verified": 373,
                "generator_additions_verified": 373,
                "ordered_image_sha256": digest,
            }
        )

    out_nodes = []
    for j, node in sorted(nodes.items()):
        orbit = [j, k.pow(j, 7), k.pow(j, 49)]
        canonical_orbit = min(tuple(orbit[i:] + orbit[:i]) for i in range(3))
        stable, mixed = counts_for_module(7, 3, node["frobenius_mod_3_representative"])
        out_nodes.append(
            {
                **node,
                "rational_3_kernel_candidates": stable,
                "mixed_3_kernel_candidates": mixed,
                "j_frobenius_orbit": canonical_orbit,
                "endomorphism_discriminant_inferred": -59 if mixed else -531,
            }
        )

    # A rejected direct screen is repaired by a single actual isogeny. The
    # surface models are Galois/scaling-related to the earlier mapped curve.
    prior_path = (
        Path(__file__).resolve().parent / "results/quotient-divisor-verification.json"
    )
    prior = json.loads(prior_path.read_text())
    anchor = prior["sources"][0]
    if (anchor["source_a"], anchor["source_b"], anchor["source_points_lifted"]) != (
        248,
        42,
        373,
    ):
        raise ArithmeticError("Unexpected prior divisor-map certificate")
    target_scalar = anchor["source_generator_image_scalar_in_target_generator"]
    routes = []
    for node in out_nodes:
        if node["mixed_3_kernel_candidates"]:
            continue
        edge = next(e for e in edges if e["source_j"] == node["j"])
        target = nodes[edge["target_j"]]
        conjugation = next(
            i for i in range(3) if k.pow(first["j"], 7**i) == target["j"]
        )
        u = scaling(
            k,
            (k.pow(first["a"], 7**conjugation), k.pow(first["b"], 7**conjugation)),
            (target["a"], target["b"]),
        )

        # Target -> old source: undo scaling, then inverse coefficient Frobenius.
        def back_to_original(point, u=u, conjugation=conjugation):
            if point is None:
                return None
            return (
                k.pow(k.div(point[0], k.pow(u, 2)), 7 ** ((3 - conjugation) % 3)),
                k.pow(k.div(point[1], k.pow(u, 3)), 7 ** ((3 - conjugation) % 3)),
            )

        g = tuple(edge["source_generator"])
        h = back_to_original(tuple(edge["image_generator"]))
        original_generator, current = (0, 176), None
        # Scalar labels are witnesses for comparison with PR6, never an
        # inversion algorithm or part of the source-to-source map.
        label = None
        for n in range(373):
            if current == h:
                label = n
            current = ec_add(k, 0, first["a"], current, original_generator)
        if label in (None, 0):
            raise ArithmeticError("Transported generator is missing")
        current = None
        for n in range(373):
            image = evaluate_isogeny(
                k,
                node["a"],
                node["b"],
                edge["kernel_x"],
                current,
                edge["codomain_scaling"],
            )
            if back_to_original(image) != ec_mul(k, 0, first["a"], h, n):
                raise ArithmeticError("Composed route failed")
            current = ec_add(k, 0, node["a"], current, g)
        routes.append(
            {
                "source_j": node["j"],
                "surface_j": target["j"],
                "isogeny_kernel_x": edge["kernel_x"],
                "surface_conjugation_power": conjugation,
                "surface_scale": u,
                "source_generator": g,
                "original_source_image": h,
                "image_scalar_relative_to_PR6_generator": label,
                "inherited_target_scalar_relative_to_PR6_J_generator": target_scalar
                * label
                % 373,
                "composed_point_checks": 373,
                "prior_divisor_report": "quotient-divisor-verification.json",
                "new_Jacobian_inversion_performed": False,
            }
        )

    # D_pi=-531=-59*3^2. Scalar -I means (pi+1)/3 is an integral
    # endomorphism; its discriminant is -59, already fundamental.
    # The only other possible order has discriminant -531.
    return {
        "scope": {
            "toy_only": True,
            "ell": 3,
            "explicit_isogeny_transport": True,
            "mixed_kernel_screen_is_construction_specific": True,
            "new_Jacobian_maps_evaluated": False,
            "inherited_map_uses_prior_enumerative_inversion": True,
            "practical_speedup_established": False,
            "novel_vulnerability_established": False,
        },
        "field": {
            "p": 7,
            "degree": 3,
            "modulus_ascending": k.modulus,
            "encoding": "a0+7*a1+49*a2",
        },
        "prior_divisor_report_sha256": sha256(prior_path.read_bytes()).hexdigest(),
        "endomorphism_inference": {
            "ordinary": True,
            "trace": -29,
            "frobenius_discriminant": -531,
            "fundamental_discriminant": -59,
            "frobenius_order_conductor": 3,
            "scalar_endomorphism_trace": -9,
            "scalar_endomorphism_norm": 35,
        },
        "summary": {
            "j_vertices": len(nodes),
            "directed_edges": len(edges),
            "vertices_by_mixed_kernel_count": dict(
                Counter(row["mixed_3_kernel_candidates"] for row in out_nodes)
            ),
            "j_frobenius_orbits": sorted(
                {row["j_frobenius_orbit"] for row in out_nodes}
            ),
            "point_images_verified": 373 * len(edges),
            "generator_additions_verified": 373 * len(edges),
        },
        "nodes": out_nodes,
        "edges": edges,
        "routes_to_prior_descent_example": routes,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
