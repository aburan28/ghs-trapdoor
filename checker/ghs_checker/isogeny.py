"""Explicit K-rational cyclic 3-isogenies in characteristic two.

Kernel x-coordinates solve u^4+u^3+b=0. The kernel is Frobenius-stable
when u is in K, even if the individual nonzero kernel points are not.
"""

from collections import deque
from dataclasses import dataclass
from .analysis import analyze_curve
from .curve import BinaryCurve


def _trim(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def _add(a, b):
    result = [0] * max(len(a), len(b))
    for i, value in enumerate(a):
        result[i] ^= value
    for i, value in enumerate(b):
        result[i] ^= value
    return _trim(result)


def _mod(a, b, field):
    a = a.copy()
    if not b:
        raise ZeroDivisionError("Zero polynomial divisor")
    inv = field.inv(b[-1])
    while a and len(a) >= len(b):
        shift, factor = len(a) - len(b), field.mul(a[-1], inv)
        for i, value in enumerate(b):
            if value:
                a[i + shift] ^= field.mul(factor, value)
        _trim(a)
    return a


def _quotient(a, b, field):
    a, quotient = a.copy(), [0] * max(0, len(a) - len(b) + 1)
    inv = field.inv(b[-1])
    while a and len(a) >= len(b):
        shift, factor = len(a) - len(b), field.mul(a[-1], inv)
        quotient[shift] = factor
        for i, value in enumerate(b):
            if value:
                a[i + shift] ^= field.mul(factor, value)
        _trim(a)
    if a:
        raise ArithmeticError("Inexact polynomial division")
    return _trim(quotient)


def _gcd(a, b, field):
    while b:
        a, b = b, _mod(a, b, field)
    if a:
        inv = field.inv(a[-1])
        a = [field.mul(c, inv) for c in a]
    return a


def _square_mod(a, modulus, field):
    result = [0] * (2 * len(a) - 1) if a else []
    for i, value in enumerate(a):
        result[2 * i] = field.square(value)
    return _mod(_trim(result), modulus, field)


def _eval(poly, x, field):
    value = 0
    for coefficient in reversed(poly):
        value = field.mul(value, x) ^ coefficient
    return value


def linear_roots(poly, field, exhaustive=None):
    """All roots in K, not just roots whose associated kernel points are K-rational."""
    poly = _trim(poly.copy())
    if exhaustive is None:
        exhaustive = field.degree <= 10
    if exhaustive:
        if field.degree > 16:
            raise ValueError("Exhaustive root search is limited to degree <=16")
        return [x for x in range(field.size) if _eval(poly, x, field) == 0]
    x = [0, 1]
    xq = x
    for _ in range(field.degree):
        xq = _square_mod(xq, poly, field)
    rational = _gcd(poly, _add(xq, x), field)
    roots = []

    def split(factor):
        if len(factor) <= 1:
            return
        if len(factor) == 2:
            roots.append(field.div(factor[0], factor[1]))
            return
        # Trace pairings with a polynomial basis separate any two distinct
        # elements of K. This deterministic loop avoids unbounded random retries.
        for bit in range(field.degree):
            power, trace = [0, 1 << bit], []
            power = _mod(power, factor, field)
            for _ in range(field.degree):
                trace = _add(trace, power)
                power = _square_mod(power, factor, field)
            part = _gcd(factor, trace, field)
            if 1 < len(part) < len(factor):
                split(part)
                split(_quotient(factor, part, field))
                return
        raise ArithmeticError("Failed to split the rational-root factor")

    split(rational)
    return sorted(roots)


@dataclass(frozen=True)
class ThreeIsogeny:
    source: BinaryCurve
    kernel_x: int

    def __post_init__(self):
        f, u = self.source.field, self.kernel_x
        f.element(u)
        if f.square(f.square(u)) ^ f.mul(f.square(u), u) ^ self.source.b:
            raise ValueError("Kernel x-coordinate does not satisfy the 3-division polynomial")

    @property
    def target(self):
        f, u = self.source.field, self.kernel_x
        return BinaryCurve(f, self.source.a, self.source.b ^ u ^ f.square(u))

    def apply(self, point):
        """Vélu map followed by the y-translation normalizing a4 to zero."""
        if not self.source.contains(point):
            raise ValueError("Point is not on the source curve")
        if point is None:
            return None
        f, u = self.source.field, self.kernel_x
        x, y = point
        if x == u:
            return None
        inv = f.inv(x ^ u)
        inv2, x2 = f.square(inv), f.square(x)
        x3 = f.mul(x2, x)
        xx = x ^ f.mul(f.mul(u, x), inv2)
        numerator = x3 ^ f.mul(u ^ 1, x2)
        yy = y ^ f.mul(f.mul(u, y), inv2) ^ f.mul(f.mul(u, numerator), f.mul(inv2, inv))
        result = xx, yy
        if not self.target.contains(result):
            raise ArithmeticError("Isogeny image failed the target curve equation")
        return result

    def describe(self):
        return {"degree": 3, "kernel_x": hex(self.kernel_x), "target_b": hex(self.target.b),
                "kernel_polynomial": "X+kernel_x", "rational_over_input_field": True,
                "individual_kernel_points_need_not_be_rational": True,
                "map": "ghs_checker.isogeny.ThreeIsogeny.apply"}


def three_isogenies(curve):
    return [ThreeIsogeny(curve, u) for u in linear_roots([curve.b, 0, 0, 1, 1], curve.field)]


def scan_isogenies(curve, subgroup_order=None, subfields=None, max_genus=32,
                   max_curves=8, max_depth=3, metadata=None):
    if not isinstance(max_curves, int) or isinstance(max_curves, bool) or not 1 <= max_curves <= 10000:
        raise ValueError("max_curves must be between 1 and 10000")
    if not isinstance(max_depth, int) or isinstance(max_depth, bool) or not 0 <= max_depth <= 100:
        raise ValueError("max_depth must be between 0 and 100")
    if subgroup_order is not None and subgroup_order % 3 == 0:
        raise ValueError("This scanner requires a subgroup order coprime to 3; other cases need point-specific injectivity checks")
    root = analyze_curve(curve, subgroup_order, subfields, max_genus, metadata)
    root.update(search={"degrees": [3], "max_curves": max_curves, "max_depth": max_depth,
                        "bounded_search_only": True,
                        "other_isogeny_degrees_searched": False,
                        "subgroup_injectivity_of_isogeny_paths": "conditional_on_claimed_order_coprime_to_3" if subgroup_order else "not_checked"})
    nodes, edges = [], []
    seen = {curve.b: 0}
    queue = deque([(curve, 0, [])])
    budget_hit = False
    depth_hit = False
    while queue:
        current, depth, path = queue.popleft()
        index = seen[current.b]
        report = analyze_curve(current, subgroup_order, subfields, max_genus)
        nodes.append({"index": index, "depth": depth, "curve": report["curve"],
                      "path_kernel_x": path, "descents": report["descents"]})
        if depth >= max_depth:
            depth_hit = True
            continue
        for phi in three_isogenies(current):
            target = phi.target
            if target.b not in seen:
                if len(seen) >= max_curves:
                    budget_hit = True
                    continue
                seen[target.b] = len(seen)
                entry = target, depth + 1, path + [hex(phi.kernel_x)]
                queue.append(entry)
            edges.append({"source": index, "target": seen[target.b], **phi.describe()})
    root["nodes"] = sorted(nodes, key=lambda node: node["index"])
    root["edges"] = edges
    root["search"].update(curves_checked=len(nodes), edges_recorded=len(edges),
                          curve_budget_reached=budget_hit, depth_limit_reached=depth_hit,
                          scope="Only the explored K-rational cyclic 3-isogeny component within the stated bounds")
    root["candidate_nodes"] = [node["index"] for node in root["nodes"] if any(
        d["screening"] == "small_genus_candidate" for d in node["descents"])]
    root["interpretation"] = "Candidates require an actual descent map and cost analysis; absence of candidates is not a security certificate"
    return root
