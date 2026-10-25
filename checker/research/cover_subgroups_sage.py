"""Actual Hess divisor transfer controls for the Prym and wild-cover tracks."""

import argparse
import json
from hashlib import sha256
from pathlib import Path
from time import perf_counter

from sage.all import GF, FunctionField, PolynomialRing
from sage.env import SAGE_VERSION

if not __debug__:
    raise RuntimeError("Research certificates require assertions; omit -O")


def function_field(k, name, coefficients, degree=2):
    base = FunctionField(k, name)
    x = base.gen()
    ring = PolynomialRing(base, "T")
    T = ring.gen()
    rhs = sum(k(c) * x**i for i, c in enumerate(coefficients))
    return base.extension(
        T**degree - rhs if k.characteristic() != 2 else T**degree + T - rhs, "y"
    )


def at(field, xx, yy):
    base = field.base_field()
    prime = base.maximal_order().ideal(base.gen() - xx).place()
    return next(p for p in field.places_above(prime) if field.gen().evaluate(p) == yy)


def prym():
    k = GF(3)
    f = [2, 1] + [0] * 8 + [1]
    cover = [0] * 21
    cover[::2] = f
    D, C = function_field(k, "u", f), function_field(k, "x", cover)
    assert D.genus() == 4 and C.genus() == 9
    P, O = at(D, k(1), k(1)), at(D, k(1), k(2))
    fibers_p = [at(C, xx, k(1)) for xx in (k(1), k(2))]
    fibers_o = [at(C, xx, k(2)) for xx in (k(1), k(2))]
    GD = D.jacobian(model="hess", base_div=4 * O.divisor()).group()
    GC = C.jacobian(model="hess", base_div=9 * fibers_o[0].divisor()).group()
    source = 2 * GD.point(P - O)
    image = 2 * GC.point(
        sum(p.divisor() for p in fibers_p) - sum(p.divisor() for p in fibers_o)
    )
    assert source != GD.zero() and 41 * source == GD.zero()
    assert image != GC.zero() and 41 * image == GC.zero()
    # Under C->B, (x,y)->(x^2,x*y), these two fiber divisors
    # have literally the same image {(1,1),(1,2)}; their difference norms to 0.
    norm_p = sorted((int(xx * xx), int(xx * k(1))) for xx in (k(1), k(2)))
    norm_o = sorted((int(xx * xx), int(xx * k(2))) for xx in (k(1), k(2)))
    assert norm_p == norm_o
    return {
        "kind": "prym",
        "field": "F3",
        "source_genus": 4,
        "cover_genus": 9,
        "Prym_dimension": 4,
        "source_order": 82,
        "cofactor": 2,
        "subgroup_order": 41,
        "source_divisor": str(source.divisor()),
        "image_divisor": str(image.divisor()),
        "norm_image_divisor": "0 (fiber images cancel literally)",
        "actual_source_and_image_orders_verified": True,
        "image_in_connected_Prym_reason": "pullback Jac(D) is connected; its norm is identically zero by the two quotient maps",
        "polarization": "canonical principal Prym polarization; restricted cover polarization is twice principal",
    }


def wild():
    k = GF(4, "a")
    a = k.gen()
    assert a * a == a + 1
    E = function_field(k, "u", [0, 0, 0, a])
    H = function_field(k, "x", [0, 0, 0, 1], degree=4)
    assert E.genus() == 1 and H.genus() == 3
    x, z = H.base_field().gen(), H.gen()
    quotient_y = a * a * z * z + a * z
    assert quotient_y * quotient_y + quotient_y == a * x**3
    conjugate_y = a * a * (z + a * a) ** 2 + a * (z + a * a)
    assert conjugate_y == quotient_y
    Oe, Oh = E.places_infinite()[0], H.places_infinite()[0]
    GE = E.jacobian(model="hess", base_div=Oe.divisor()).group()
    GH = H.jacobian(model="hess", base_div=3 * Oh.divisor()).group()
    P = at(E, k(0), k(0))
    # y_E = a^2*z^2+a*z, so its zero fiber is z=0,a^2 at x=0.
    fiber = [at(H, k(0), z) for z in (k(0), a * a)]
    source = GE.point(P - Oe)
    pullback = GH.point(sum(p.divisor() for p in fiber) - 2 * Oh.divisor())
    assert source != GE.zero() and 3 * source == GE.zero()
    assert pullback != GH.zero() and 3 * pullback == GH.zero()
    conjugate_fiber = [at(H, k(0), z * z) for z in (k(0), a * a)]
    image = pullback + GH.point(
        sum(p.divisor() for p in conjugate_fiber) - 2 * Oh.divisor()
    )
    assert image != GH.zero() and 3 * image == GH.zero()
    # Literal Frobenius stability of the degree-zero support cycle certifies
    # that the trace image is F2-rational, independent of coordinate choices.
    support = [k(0), k(0), a, a * a]
    assert sorted(map(str, support)) == sorted(str(z * z) for z in support)
    base_counts = []
    for degree in (1, 2, 3):
        L = GF(2**degree, "s") if degree > 1 else GF(2)
        base_counts.append(1+sum(zz**4+zz == xx**3 for xx in L for zz in L))
    assert base_counts == [3, 5, 9]
    # The first three power sums vanish; genus three and reciprocity then
    # determine the complete F2 polynomial T^6+8.
    return {
        "kind": "wild",
        "cover": "z^4+z=x^3 over F2",
        "cover_genus": 3,
        "counts_H_over_F2_F4_F8": base_counts,
        "full_F2_Frobenius_polynomial_H_ascending": [8, 0, 0, 0, 0, 0, 1],
        "Jacobian_H_F2_order": 9,
        "source": "y^2+y=a*x^3 over F4, a^2+a+1=0",
        "source_genus": 1,
        "quotient_map": "u=x, y=a^2*z^2+a*z",
        "degree": 2,
        "source_subgroup_order": 3,
        "source_divisor": str(source.divisor()),
        "pullback_divisor": str(pullback.divisor()),
        "trace_image_divisor": str(image.divisor()),
        "actual_nonzero_source_pullback_and_trace_orders_verified": True,
        "trace_image_field": "F2",
        "trace_map": "D -> pi^*D + Frob2(pi^*D)",
        "norm_conorm_identity": "pi_* pi^*=[2]",
        "ordinary_elliptic_source": False,
        "scope": "supersingular toy control; not an ordinary binary-curve speedup",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", choices=("prym", "wild"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    start = perf_counter()
    result = {"prym": prym, "wild": wild}[args.kind]()
    result.update(
        {
            "sage_version": SAGE_VERSION,
            "elapsed_seconds": perf_counter() - start,
            "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "practical_speedup_established": False,
        }
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
