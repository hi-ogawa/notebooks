# x^e (mod n)
def pow_mod(x, e, n):
    y = 1
    while e > 0:
        if e & 1:
            y = (y * x) % n
        x = (x * x) % n
        e = e >> 1
    return y


# x^-1 \in Zp
def inv_mod(x, p):
    return pow_mod(x, p - 2, p)


def factorize(n):
    from math import isqrt

    res = []
    for x in range(2, isqrt(n) + 1):
        while n % x == 0:
            n //= x
            res.append(x)
    if n > 1:
        res.append(n)
    return res


#
# Prime testing
#
class MillerRabin:
    def __init__(self, seed="miller-rabin", num_trials=128):
        import random

        self.rng = random.Random(seed)
        self.num_trials = num_trials

    # x^(n - 1) != 1  =>  n is composite
    def test(self, x, k, q, n):
        y = pow_mod(x, q, n)
        if y == 1:
            return False
        for _ in range(k):
            if y == n - 1:
                return False
            y = (y * y) % n
        return True

    def is_composite(self, n):
        from math import gcd

        assert n >= 0
        if n == 1 or n % 2 == 0:
            return False
        if n == 2:
            return True

        # n = 2^k q + 1
        q = n - 1
        k = 0
        while q % 2 == 0:
            q //= 2
            k += 1

        for _ in range(self.num_trials):
            x = self.rng.randrange(1, n)
            if gcd(x, n) != 1:
                return True
            if self.test(x, k, q, n):
                return True
        return False

    def is_prime(self, n):
        return not self.is_composite(n)


MILLER_RABIN = MillerRabin()

#
# Sqrt mod
#
class TonelliShanks:
    def __init__(self, seed="tonelli-shanks"):
        import random

        self.rng = random.Random(seed)  # For finding non quadratic residue

    # x^2 = y ∈ Z(p)
    def solve(self, y, p):
        y = y % p
        if y == 0 or p == 2:
            return y

        # p = 2^k q + 1 (odd prime)
        q = p - 1
        k = 0
        while q % 2 == 0:
            q //= 2
            k += 1

        # Euler's criterion
        if pow_mod(y, (p - 1) // 2, p) != 1:
            return None

        # Find non quadratic residue
        # z^((p - 1)/2) = -1
        while True:
            z = self.rng.randrange(1, p)
            if pow_mod(z, (p - 1) // 2, p) != 1:
                break

        # Tonelli Shanks algorithm

        # y^(q + 1) = y^q y = a y
        #           = (y^(q + 1)/2)^2 = b^2
        a = pow_mod(y, q, p)
        b = pow_mod(y, (q + 1) // 2, p)

        # y = (a^-1) b^2 = ((a^-1/2) b)^2
        # (accumulate the factor "a^-1/2" during the loop below)
        res = b

        # a^(2^(k - 1)) = y^(q 2^(k - 1)) = y^((n - 1) / 2) = 1
        while a != 1:
            # Find t = min{i | a^(2^i) = 1} thus
            #   a^(2^t)     = 1
            #   a^(2^(t-1)) = -1
            t = None
            d = a
            for i in range(1, k):
                d = (d * d) % p
                if d == 1:
                    t = i
                    break
            assert t is not None

            # Update "a" which leads to smaller "t" on the next loop
            # 1 = (-1)^2
            #   = a^(2^(t-1)) z^((p-1)/2)
            #   = a^(2^(t-1)) z^(q 2^(k-1))
            #   = a^(2^(t-1)) (z^(q 2^(k-t-1))^2)^(2^(t-1))
            #   = (a w^2)^(2^(t-1))
            w = pow_mod(z, q * (2 ** (k - t - 1)), p)
            a = a * w % p * w % p
            res = res * w % p

        # Pick smaller one
        return min(res, p - res)


TONELLI_SHANKS = TonelliShanks()


class PointBase:
    """
    # y^2 = x^3 - x in Z61
    >>> Point = make_point_class(-1, 0, 61)
    >>> p = Point.lift(4); p
    (4, 11)
    >>> Point.zero()
    oo
    >>> order = 1 # infinity
    >>> for x in range(61):
    ...     y = Point.x_to_y(x)
    ...     if y is None:
    ...         continue
    ...     order += 1 + (y != 0)
    >>> order
    72
    """

    a = None
    b = None
    m = None

    #
    # Modulo integer routines
    #
    @classmethod
    def mod(cls, x):
        return x % cls.m

    @classmethod
    def mod_pow(cls, x, e):
        y = 1
        while e > 0:
            if e & 1:
                y = cls.mod(y * x)
            x = cls.mod(x * x)
            e = e >> 1
        return y

    @classmethod
    def mod_inv(cls, x):
        return cls.mod_pow(x, cls.m - 2)

    @classmethod
    def mod_div(cls, x, y):
        return cls.mod(x * cls.mod_inv(y))

    #
    # Elliptic curve
    #
    @classmethod
    def rhs(cls, x):
        return cls.mod(x ** 3 + cls.a * x + cls.b)

    @classmethod
    def x_to_y(cls, x):
        return TONELLI_SHANKS.solve(cls.rhs(x), cls.m)

    @classmethod
    def lift(cls, x):
        y = TONELLI_SHANKS.solve(cls.rhs(x), cls.m)
        assert y is not None
        return cls(x, y)

    @classmethod
    def zero(cls):
        u = cls(0, 0)
        u.is_zero = True
        return u

    #
    # Methods
    #
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.is_zero = False

    def copy(self):
        if self.is_zero:
            return type(self).zero()
        return type(self)(*self)

    def check(self):
        x, y = self
        return self.mod(y ** 2) == self.rhs(x)

    def __repr__(self):
        if self.is_zero:
            return "oo"
        return f"({self.x}, {self.y})"

    def __eq__(self, q):
        p = self
        return type(p) == type(q) and self.is_zero == q.is_zero and tuple(p) == tuple(q)

    #
    # Destructuring by "x, y = point"
    #
    def __len__(self):
        return 2

    def __iter__(self):
        yield self.x
        yield self.y

    #
    # Group structure
    #
    def add_general(self, q):
        p = self
        assert not (p.is_zero or q.is_zero or p.x == q.x)
        x1, y1 = p
        x2, y2 = q
        s = self.mod_div(y2 - y1, x2 - x1)
        x = self.mod(s ** 2 - (x1 + x2))
        y = self.mod(-s * (x - x1) - y1)
        return type(self)(x, y)

    def double_general(self):
        assert not (self.is_zero or self.y == 0)
        x1, y1 = self
        s = self.mod_div(3 * x1 ** 2 + self.a, 2 * y1)
        x = self.mod(s ** 2 - 2 * x1)
        y = self.mod(-s * (x - x1) - y1)
        return type(self)(x, y)

    def double(self):
        if self.is_zero or self.y == 0:
            return type(self).zero()
        return self.double_general()

    def mul(self, c):
        assert isinstance(c, int)
        p = self
        q = self.zero()
        if c < 0:
            c = -c
            p = -p
        while c > 0:
            if c & 1:
                q = q + p
            p = p.double()
            c = c >> 1
        return q

    def __add__(self, q):
        """
        >>> Point = make_point_class(-1, 0, 61)
        >>> z = Point.zero()
        >>> p = Point.lift(4); p
        (4, 11)
        >>> q = Point.lift(6); q
        (6, 24)
        >>> p + q
        (17, 57)
        >>> q + p
        (17, 57)
        >>> p + z
        (4, 11)
        >>> z + p
        (4, 11)
        >>> q + q
        (46, 19)
        """
        p = self
        if p.is_zero:
            return q
        if q.is_zero:
            return p
        if p == q:
            return p.double()
        if p.x == q.x:
            return type(self).zero()
        return p.add_general(q)

    def __neg__(self):
        """
        >>> Point = make_point_class(-1, 0, 61)
        >>> p = Point.lift(4)
        >>> -p
        (4, 50)
        >>> p - p
        oo
        """
        p = self.copy()
        p.y = self.m - p.y
        return p

    def __sub__(self, q):
        return self + (-q)

    def __mul__(self, c):
        """
        >>> Point = make_point_class(-1, 0, 61)
        >>> p = Point.lift(6)
        >>> p
        (6, 24)
        >>> 2 * p
        (46, 19)
        >>> 3 * p
        (50, 12)
        >>> 4 * p
        (57, 1)
        >>> for i in range(1, 100):
        ...     p = i * Point.lift(6)
        ...     print(i, p)
        ...     if p == Point.zero():
        ...         break
        1 (6, 24)
        2 (46, 19)
        3 (50, 12)
        4 (57, 1)
        5 (10, 40)
        6 (0, 0)
        7 (10, 21)
        8 (57, 60)
        9 (50, 49)
        10 (46, 42)
        11 (6, 37)
        12 oo
        >>> for i in range(1, 100):
        ...     p = i * Point.lift(4)
        ...     print(i, p)
        ...     if p == Point.zero():
        ...         break
        1 (4, 11)
        2 (4, 50)
        3 oo
        """
        if not isinstance(c, int):
            return NotImplemented
        return self.mul(c)

    def __rmul__(self, c):
        return self * c


def make_point_class(a, b, m):
    assert MILLER_RABIN.is_prime(m)
    return type(f"Point[a = {a}, b = {b}, p = {m}]", (PointBase,), dict(a=a, b=b, m=m))


def demo1():
    print("\n--- demo1 ---\n")
    miller_rabin = MillerRabin()
    for e in range(2, 32):
        n = 2 ** e - 1
        print(e, n, factorize(n), miller_rabin.is_prime(n))


def demo2():
    print("\n--- demo2 ---\n")
    a = 0
    b = 7
    modulo = 2 ** 256 - 2 ** 32 - 977
    assert MILLER_RABIN.is_prime(modulo)

    Point = make_point_class(a, b, modulo)
    print(Point)

    gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
    g = Point.lift(gx)
    assert g.check()

    print(g)
    print(2 * g)
    print(-g)


if __name__ == "__main__":
    demo1()
    demo2()
