"""Small odd finite fields and exact arithmetic for bounded research experiments.

An integer sum(a_i*p**i) encodes sum(a_i*z**i). Moduli and polynomial
coefficients are in ascending order. This intentionally supports only degrees
1..3 and small primes; it is not a cryptographic finite-field implementation.
"""

from itertools import product

from .genus2_audit import prime64


class SmallField:
    def __init__(self, p, degree=1, modulus=None):
        if type(p) is not int or not 5 <= p <= 31 or not prime64(p):
            raise ValueError('requires a prime between 5 and 31')
        if type(degree) is not int or degree not in (1, 2, 3):
            raise ValueError('requires extension degree 1, 2 or 3')
        self.p, self.degree, self.size = p, degree, p**degree
        if modulus is None:
            modulus = next(list(c)+[1] for c in product(range(p), repeat=degree)
                           if degree == 1 or not any(
                               (sum(c[i]*x**i for i in range(degree))+x**degree) % p == 0
                               for x in range(p)))
        if (len(modulus) != degree+1 or modulus[-1] != 1 or
                any(type(c) is not int or not 0 <= c < p for c in modulus)):
            raise ValueError('modulus must be monic, with coefficients in [0,p)')
        # In degrees 2 and 3, having no root is equivalent to irreducibility.
        if degree > 1 and any(sum(c*x**i for i, c in enumerate(modulus)) % p == 0
                              for x in range(p)):
            raise ValueError('modulus is reducible')
        self.modulus = tuple(modulus)
        self._digits = [tuple((x//p**i) % p for i in range(degree))
                        for x in range(self.size)]
        self._add = self._mul = None

    def elt(self, *coefficients):
        if len(coefficients) > self.degree:
            raise ValueError('too many coefficients')
        return sum((c % self.p)*self.p**i for i, c in enumerate(coefficients))

    def digits(self, x):
        return self._digits[x]

    def add(self, a, b):
        if self._add is not None:
            return self._add[a][b]
        return self.elt(*(x+y for x, y in zip(self._digits[a], self._digits[b])))

    def neg(self, a):
        return self.elt(*(-x for x in self._digits[a]))

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        if self._mul is not None:
            return self._mul[a][b]
        p, d = self.p, self.degree
        u, v = self._digits[a], self._digits[b]
        c = [sum(u[j]*v[i-j] for j in range(d) if 0 <= i-j < d)
             for i in range(2*d-1)]
        for i in range(len(c)-1, d-1, -1):
            for j in range(d):
                c[i-d+j] -= c[i]*self.modulus[j]
        return self.elt(*c[:d])

    def pow(self, a, n):
        if n < 0:
            return self.pow(self.inv(a), -n)
        out = 1
        while n:
            if n & 1:
                out = self.mul(out, a)
            a, n = self.mul(a, a), n >> 1
        return out

    def inv(self, a):
        if not a:
            raise ZeroDivisionError
        return self.pow(a, self.size-2)

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def tables(self):
        """Optional tiny-field acceleration; never allocate quadratic large tables."""
        if self.size > 343:
            raise ValueError('tables are limited to fields of size <=343')
        if self._add is None:
            additions = [[self.add(a, b) for b in range(self.size)]
                         for a in range(self.size)]
            products = [[self.mul(a, b) for b in range(self.size)]
                        for a in range(self.size)]
            self._add, self._mul = additions, products
        return self._add, self._mul


def factors(n):
    result, d = [], 2
    while d*d <= n:
        exponent = 0
        while n % d == 0:
            n, exponent = n//d, exponent+1
        if exponent:
            result.append([d, exponent])
        d += 1
    if n > 1:
        result.append([n, 1])
    return result


def evaluate(k, coefficients, x):
    out = 0
    for c in reversed(coefficients):
        out = k.add(k.mul(out, x), c)
    return out


def points(k, rhs):
    """Affine points; coefficients already encoded in k."""
    roots = {}
    for y in range(k.size):
        roots.setdefault(k.mul(y, y), []).append(y)
    for x in range(k.size):
        for y in roots.get(evaluate(k, rhs, x), ()):
            yield x, y


def curve_count(k, rhs):
    """Smooth hyperelliptic projective count, including both charts at infinity."""
    if len(rhs) < 2 or not rhs[-1]:
        raise ValueError('positive degree polynomial with nonzero leading coefficient required')
    squares = {k.mul(y, y) for y in range(k.size)}
    n = sum(1 if (value := evaluate(k, rhs, x)) == 0 else
            2 if value in squares else 0 for x in range(k.size))
    return n + (1 if (len(rhs)-1) % 2 else 2 if rhs[-1] in squares else 0)


def ec_add(k, a2, a4, left, right):
    """Chord law on y^2=x^3+a2*x^2+a4*x+a6, O=None; smooth points only."""
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if k.add(y, v) == 0:
            return None
        numerator = k.add(k.add(k.mul(3, k.mul(x, x)),
                                k.mul(k.mul(2, a2), x)), a4)
        slope = k.div(numerator, k.mul(2, y))
    else:
        slope = k.div(k.sub(v, y), k.sub(u, x))
    z = k.sub(k.sub(k.sub(k.mul(slope, slope), a2), x), u)
    return z, k.sub(k.mul(slope, k.sub(x, z)), y)


def ec_mul(k, a2, a4, point, n):
    if n < 0:
        return ec_mul(k, a2, a4, None if point is None else
                      (point[0], k.neg(point[1])), -n)
    out = None
    while n:
        if n & 1:
            out = ec_add(k, a2, a4, out, point)
        point, n = ec_add(k, a2, a4, point, point), n >> 1
    return out


def frobenius_polynomial(p, counts):
    """P(T)=det(T-Frob), ascending coefficients, from g counts for genus g."""
    g = len(counts)
    sums = [p**n+1-counts[n-1] for n in range(1, g+1)]
    descending = [1]
    for n in range(1, g+1):
        value = -sum(descending[n-j]*sums[j-1] for j in range(1, n+1))
        if value % n:
            raise ArithmeticError('Newton coefficient is not integral')
        descending.append(value//n)
    out = [0]*(2*g+1)
    for i, c in enumerate(descending):
        out[2*g-i], out[i] = c, c*p**(g-i)
    return out


def integer_poly_mul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out
