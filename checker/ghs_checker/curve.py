"""Ordinary binary elliptic curves y^2+xy=x^3+a*x^2+b."""

from dataclasses import dataclass
from .field import BinaryField


@dataclass(frozen=True)
class BinaryCurve:
    field: BinaryField
    a: int
    b: int

    def __post_init__(self):
        self.field.element(self.a)
        self.field.element(self.b)
        if not self.b:
            raise ValueError("b=0 gives a singular curve")

    def contains(self, point):
        if point is None:
            return True
        x, y = point
        if not (0 <= x < self.field.size and 0 <= y < self.field.size):
            return False
        f = self.field
        return f.square(y) ^ f.mul(x, y) == f.mul(f.square(x), x) ^ f.mul(self.a, f.square(x)) ^ self.b

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        f = self.field
        x, y = p
        xx, yy = q
        if x == xx:
            if y != yy or x == 0:
                return None
            slope = x ^ f.div(y, x)
            xxx = f.square(slope) ^ slope ^ self.a
            yyy = f.square(x) ^ f.mul(slope ^ 1, xxx)
        else:
            slope = f.div(y ^ yy, x ^ xx)
            xxx = f.square(slope) ^ slope ^ x ^ xx ^ self.a
            yyy = f.mul(slope, x ^ xxx) ^ xxx ^ y
        return xxx, yyy

    def multiply(self, scalar, point):
        if scalar < 0:
            scalar, point = -scalar, None if point is None else (point[0], point[0] ^ point[1])
        result = None
        while scalar:
            if scalar & 1:
                result = self.add(result, point)
            point = self.add(point, point)
            scalar >>= 1
        return result

    def count_points(self):
        """Exhaustive x/trace point count, deliberately restricted to small fields."""
        if self.field.degree > 16:
            raise ValueError("Built-in exact point counting is limited to degree <=16")
        f = self.field
        count = 2  # infinity, and the unique point with x=0
        for x in range(1, f.size):
            value = x ^ self.a ^ f.div(self.b, f.square(x))
            if f.absolute_trace(value) == 0:
                count += 2
        return count

    def describe(self):
        return {"model": "y^2+x*y=x^3+a*x^2+b", "a": hex(self.a), "b": hex(self.b),
                "j": hex(self.field.inv(self.b)), "nonsingular": True}
