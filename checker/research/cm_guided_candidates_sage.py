"""Sage 10.6: bounded CM-guided cubic-field kernel candidates.

This constructs actual curves and torsion data, not quotient Jacobians or
cryptographic vulnerabilities. Run with sage -python and --output PATH.
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    GF,
    ZZ,
    EllipticCurve,
    PolynomialRing,
    QuadraticField,
    hilbert_class_polynomial,
)
from sage.env import SAGE_VERSION

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")

CASES = ((13, -83, -211), (13, 61, -563), (19, -137, -963))
MODULI = {13: [11, 2, 0, 1], 19: [17, 4, 0, 1]}


def run():
    classes, candidates = [], []
    for p, trace, disc in CASES:
        base = GF(p)
        ring = PolynomialRing(base, "t")
        K = GF(p**3, "z", modulus=ring(MODULI[p]))
        z = K.gen()

        def encode(a, p=p):
            return sum(int(a.polynomial()[i]) * p**i for i in range(3))

        def decode(n, p=p, K=K, z=z):
            return sum(K((n // p**i) % p) * z**i for i in range(3))

        order = p**3 + 1 - trace
        assert ZZ(order).is_prime() and trace % p != 0
        assert (p**3 + trace + 1) % 9 == 0
        assert trace**2 - 4 * p**3 == 9 * disc
        H = hilbert_class_polynomial(disc)
        roots = sorted(H.change_ring(K).roots(multiplicities=False), key=encode)
        seen, selected = set(), []
        nonsquare = next(decode(n) for n in range(2, p**3) if not decode(n).is_square())
        for j in roots:
            orbit = tuple(sorted(encode(j ** (p**i)) for i in range(3)))
            if j**p == j or orbit in seen:
                continue
            seen.add(orbit)
            aa, bb = 3 * j / (K(1728) - j), 2 * j / (K(1728) - j)
            models = (
                EllipticCurve(K, [aa, bb]),
                EllipticCurve(K, [aa * nonsquare**2, bb * nonsquare**3]),
            )
            valid = [E for E in models if E.cardinality() == order]
            assert len(valid) == 1
            E = valid[0]
            division_roots = sorted(
                E.division_polynomial(3).roots(multiplicities=False), key=encode
            )
            assert len(division_roots) == 4
            assert all(
                not (x**3 + E.a4() * x + E.a6()).is_square() for x in division_roots
            )
            # Four rational x-coordinates, with all y-coordinates outside K,
            # verify scalar -I, rather than inferring it from trace alone.
            generator = None
            for n in range(p**3):
                x = decode(n)
                value = x**3 + E.a4() * x + E.a6()
                if value.is_square():
                    yy = value.sqrt()
                    yy = min((yy, -yy), key=encode)
                    generator = E(x, yy)
                    break
            assert generator is not None and not generator.is_zero()
            assert order * generator == E.zero()
            record = {
                "p": p,
                "extension_degree": 3,
                "modulus_ascending": MODULI[p],
                "a": encode(E.a4()),
                "b": encode(E.a6()),
                "j": encode(j),
                "j_frobenius_orbit_sorted": orbit,
                "group_order": order,
                "trace": trace,
                "CM_order_used_discriminant": disc,
                "endomorphism_field_discriminant": int(
                    QuadraticField(disc).discriminant()
                ),
                "endomorphism_ring_fully_determined": False,
                "three_division_roots": list(map(encode, division_roots)),
                "q3_Frobenius_mod_3_verified": [[2, 0], [0, 2]],
                "subgroup_generator": [encode(generator[0]), encode(generator[1])],
                "generator_order_verified": order,
            }
            candidates.append(record)
            selected.append(encode(j))
        assert selected
        classes.append(
            {
                "p": p,
                "trace": trace,
                "group_order": order,
                "CM_order_used_discriminant": disc,
                "hilbert_class_polynomial_ascending": list(map(int, H)),
                "class_polynomial_degree": int(H.degree()),
                "class_polynomial_roots_in_cubic_field": len(roots),
                "selected_j_orbit_representatives": selected,
            }
        )
    return {
        "scope": {
            "toy_only": True,
            "Sage_version_used": SAGE_VERSION,
            "bounded_cases": [list(case) for case in CASES],
            "all_class_polynomial_roots_in_stated_fields_examined": True,
            "quotient_targets_constructed": False,
            "source_to_Jacobian_maps_constructed": False,
            "practical_speedup_established": False,
            "novel_vulnerability_established": False,
        },
        "summary": {
            "selected_curve_models": len(candidates),
            "coefficient_Frobenius_j_orbits": len(candidates),
            "elliptic_isogeny_classes_over_stated_fields": len(classes),
            "endomorphism_fields_distinct_from_prior_minus_59": True,
        },
        "classes": classes,
        "candidates": candidates,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
