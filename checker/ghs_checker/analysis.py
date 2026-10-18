"""Classical GHS invariant/genus screening, not a complete security decision."""

from math import gcd, isqrt
from .field import integer, prime_divisors


SCOPE = {
    "construction": "classical characteristic-two GHS",
    "supported_curves": "ordinary binary y^2+x*y=x^3+a*x^2+b, b != 0",
    "descent_map_constructed": False,
    "target_subgroup_image_verified": False,
    "dlp_solved": False,
    "security_certification": False,
    "limitations": [
        "Does not classify arbitrary covers, generalized GHS, or all Weil-descent reductions.",
        "A favorable genus is a screening result, not evidence of a practical attack.",
        "No full Jacobian DLP time/memory estimate or subgroup-preserving descent map is computed.",
    ],
}


def proper_subfields(degree):
    return [s for s in range(1, degree) if degree % s == 0]


def binary_basis(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return basis


def in_span(value, basis):
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            return False
        value ^= basis[pivot]
    return True


def capped_pow(base, exponent, cap):
    result = 1
    while exponent:
        if exponent & 1:
            result = min(cap, result * base)
        exponent >>= 1
        if exponent:
            base = min(cap, base * base)
    return result


def subgroup_capacity(q, genus, order):
    if order is None:
        return {"status": "not_checked", "reason": "No subgroup order supplied"}
    # ceil((sqrt(q)+1)^2)^g is a conservative integer Hasse-Weil upper bound.
    root = isqrt(4 * q)
    base = q + 1 + root + (root * root != 4 * q)
    if genus == 1:
        base = q + 1 + root  # exact integral elliptic-curve upper bound
    impossible = capped_pow(base, genus, order) < order
    return {
        "status": "cannot_preserve_full_subgroup" if impossible else "not_excluded_by_size",
        "claimed_subgroup_order": str(order),
        "reason": "Hasse-Weil size bound excludes the claimed subgroup" if impossible else
                  "A necessary size condition only; the descent map and its kernel remain unchecked",
        "upper_bound_expression": f"{base}^{genus}",
    }


def analyze_descent(curve, subfield_degree, subgroup_order=None, max_genus=32):
    f = curve.field
    integer(subfield_degree, "Subfield degree")
    integer(max_genus, "Genus threshold")
    if not 1 <= subfield_degree < f.degree or f.degree % subfield_degree:
        raise ValueError("Choose a proper subfield degree dividing the field degree")
    if max_genus < 1:
        raise ValueError("Genus threshold must be positive")
    if subgroup_order is not None and subgroup_order < 2:
        raise ValueError("Subgroup order must be >=2")
    n = f.degree // subfield_degree
    marker = 1 << f.degree
    value, rows = curve.b, []
    for _ in range(n):
        rows.append(marker | value)
        value = f.frobenius(value, subfield_degree)
    basis = binary_basis(rows)
    m = len(basis)
    # Squaring on the second component is an invertible F_2-linear map, so using
    # b rather than sqrt(b) gives exactly the conventional augmented GHS rank.
    epsilon = int(in_span(marker, basis))
    regular = n % 2 == 1 or f.absolute_trace(curve.a) == 0 or m == n
    trace_b = f.relative_trace(curve.b, subfield_degree)
    result = {"subfield_degree": subfield_degree, "base_field": f"GF(2^{subfield_degree})",
              "extension_degree": n, "m": m,
              "rank_field": "GF(2)", "one_zero_in_augmented_span": bool(epsilon),
              "relative_trace_b": hex(trace_b), "regularity_condition_holds": regular,
              "genus": None, "genus_formula": None,
              "jacobian_log2_size_leading_term": None,
              "genus_threshold": max_genus,
              "subgroup_capacity": {"status": "not_checked"}}
    if not regular:
        result.update(screening="regularity_condition_failed",
                      reason="Even extension degree, absolute trace(a)=1, and m<n; the regular GHS genus result is not asserted")
        return result
    genus = (1 << (m - 1)) - epsilon
    capacity = subgroup_capacity(1 << subfield_degree, genus, subgroup_order)
    if capacity["status"] == "cannot_preserve_full_subgroup":
        screening = "subgroup_size_obstruction"
    elif genus <= max_genus:
        screening = "small_genus_candidate"
    else:
        screening = "above_genus_threshold"
    result.update(genus=genus, genus_formula=f"2^({m}-1)-{epsilon}",
                  jacobian_log2_size_leading_term=genus * subfield_degree,
                  subgroup_capacity=capacity, screening=screening,
                  reason="Classical GHS screen; threshold does not measure practical attack cost")
    return result


def analyze_curve(curve, subgroup_order=None, subfields=None, max_genus=32, metadata=None):
    integer(max_genus, "Genus threshold")
    if max_genus < 1:
        raise ValueError("Genus threshold must be positive")
    if subgroup_order is not None:
        integer(subgroup_order, "Subgroup order")
    if subgroup_order is not None and subgroup_order < 2:
        raise ValueError("Subgroup order must be >=2")
    chosen = proper_subfields(curve.field.degree) if subfields is None else sorted(set(subfields))
    analyses = [analyze_descent(curve, s, subgroup_order, max_genus) for s in chosen]
    return {"schema_version": 1, "scope": dict(SCOPE), "field": curve.field.describe(),
            "curve": curve.describe(), "metadata": metadata or {},
            "claimed_subgroup_order": None if subgroup_order is None else str(subgroup_order),
            "subgroup_order_verified": False,
            "subgroup_order_bits": None if subgroup_order is None else subgroup_order.bit_length(),
            "all_proper_subfields_checked": chosen == proper_subfields(curve.field.degree),
            "descents": analyses,
            "interpretation": "No proper subfield exists" if not chosen else
            "Results concern the stated classical construction only; no safe/vulnerable verdict is issued"}


def field_bounds(degree):
    integer(degree, "Field degree")
    if not 1 <= degree <= 4096:
        raise ValueError("Field degree must be between 1 and 4096")
    primes = prime_divisors(degree)
    prime = primes == [degree]
    result = {"degree": degree, "prime_degree": prime,
              "proper_subfield_degrees": proper_subfields(degree),
              "scope": "Menezes-Qu classical GHS field-wide bounds", "minimum_m_gt_one": None,
              "minimum_genus_for_m_gt_one": None,
              "interpretation": "No claim about arbitrary covers or generalized descent"}
    if prime and degree != 2 and gcd(2, degree) == 1:
        order, value = 1, 2 % degree
        while value != 1:
            value = (2 * value) % degree
            order += 1
        result.update(order_of_two_mod_degree=order, minimum_m_gt_one=order + 1,
                      minimum_genus_for_m_gt_one=(1 << order) - 1)
    return result
