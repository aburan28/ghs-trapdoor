"""Polynomial-basis arithmetic in GF(2^d), with checked field definitions."""

from dataclasses import dataclass
from functools import lru_cache
import re


def parse_element(value):
    """Read an integer or a polynomial; never evaluate input as Python."""
    if isinstance(value, bool):
        raise ValueError("Boolean values are not field elements")
    if isinstance(value, int):
        result = value
    elif isinstance(value, str):
        source = value.strip().lower().replace(" ", "").replace("_", "")
        if ("x" in source or "z" in source) and not source.startswith(("0x", "0b", "0o")):
            result = 0
            for term in source.replace("z", "x").split("+"):
                if term == "1":
                    power = 0
                elif term == "x":
                    power = 1
                elif re.fullmatch(r"x\^[0-9]+", term):
                    power = int(term[2:])
                else:
                    raise ValueError("Use polynomial terms such as x^7+x+1")
                if power > 4096:
                    raise ValueError("Polynomial exponents must be at most 4096")
                result ^= 1 << power
        else:
            try:
                result = int(source, 0) if source.startswith(("0x", "0b", "0o")) else int(source, 10)
            except ValueError as exc:
                raise ValueError("Expected a decimal, hexadecimal, binary, or polynomial value") from exc
    else:
        raise ValueError("Expected an integer or string")
    if result < 0:
        raise ValueError("Values must be nonnegative")
    return result


def integer(value, label):
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")
    return value


def polynomial_text(poly):
    return "+".join("1" if i == 0 else "x" if i == 1 else f"x^{i}"
                    for i in range(poly.bit_length() - 1, -1, -1) if (poly >> i) & 1) or "0"


def binary_remainder(a, modulus):
    if not modulus:
        raise ZeroDivisionError("Zero polynomial modulus")
    degree = modulus.bit_length() - 1
    while a.bit_length() - 1 >= degree:
        a ^= modulus << (a.bit_length() - 1 - degree)
    return a


def binary_gcd(a, b):
    while b:
        a, b = b, binary_remainder(a, b)
    return a


SQUARE_BYTES = tuple(sum(((v >> i) & 1) << (2 * i) for i in range(8)) for v in range(256))


def binary_square(a, modulus):
    result, shift = 0, 0
    while a:
        result |= SQUARE_BYTES[a & 255] << shift
        a >>= 8
        shift += 16
    return binary_remainder(result, modulus)


def prime_divisors(n):
    result, p = [], 2
    while p * p <= n:
        if n % p == 0:
            result.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result.append(n)
    return result


@lru_cache(maxsize=128)
def is_irreducible(poly):
    """Rabin's irreducibility test over F_2."""
    degree = poly.bit_length() - 1
    if degree < 1 or not poly & 1:
        return False
    x = binary_remainder(2, poly)
    targets = {degree // p for p in prime_divisors(degree)}
    power = x
    for i in range(1, degree + 1):
        power = binary_square(power, poly)
        if i in targets and binary_gcd(power ^ x, poly) != 1:
            return False
    return power == x


@dataclass(frozen=True)
class BinaryField:
    degree: int
    modulus: int

    def __post_init__(self):
        integer(self.degree, "Field degree")
        integer(self.modulus, "Field polynomial")
        if not 1 <= self.degree <= 4096:
            raise ValueError("Field degree must be between 1 and 4096")
        if self.modulus.bit_length() != self.degree + 1:
            raise ValueError("Field polynomial must have exactly the specified degree")
        if not is_irreducible(self.modulus):
            raise ValueError("Field polynomial is reducible over F_2")

    @property
    def size(self):
        return 1 << self.degree

    def element(self, value):
        value = parse_element(value)
        if value >= self.size:
            raise ValueError("Field element exceeds the polynomial-basis field representation")
        return value

    def mul(self, a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            a <<= 1
            b >>= 1
        return binary_remainder(result, self.modulus)

    def square(self, a):
        return binary_square(a, self.modulus)

    def pow(self, a, exponent):
        if exponent < 0:
            return self.pow(self.inv(a), -exponent)
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.square(a)
            exponent >>= 1
        return result

    def inv(self, a):
        if not a:
            raise ZeroDivisionError("Cannot invert zero")
        u, v, g, h = a, self.modulus, 1, 0
        while u != 1:
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v, g, h = v, u, h, g
                shift = -shift
            u ^= v << shift
            g ^= h << shift
        return binary_remainder(g, self.modulus)

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def frobenius(self, a, iterations=1):
        for _ in range(iterations % self.degree):
            a = self.square(a)
        return a

    def absolute_trace(self, a):
        result, value = 0, a
        for _ in range(self.degree):
            result ^= value
            value = self.square(value)
        if result not in (0, 1):
            raise ArithmeticError("Invalid absolute trace")
        return result

    def relative_trace(self, a, subfield_degree):
        if self.degree % subfield_degree:
            raise ValueError("Subfield degree must divide the field degree")
        result, value = 0, a
        for _ in range(self.degree // subfield_degree):
            result ^= value
            value = self.frobenius(value, subfield_degree)
        return result

    def describe(self):
        return {"characteristic": 2, "degree": self.degree,
                "modulus": polynomial_text(self.modulus), "modulus_hex": hex(self.modulus),
                "basis": "polynomial", "irreducibility_verified": True}
