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
class MillerRabin():
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

#
# Sqrt mod
#
class TonelliShanks():
    def __init__(self, seed="tonelli-shanks"):
        import random
        self.rng = random.Random(seed) # For finding non quadratic residue

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

        return res


def demo1():
    miller_rabin = MillerRabin()
    for e in range(2, 32):
        n = 2 ** e - 1
        print(e, n, factorize(n), miller_rabin.is_prime(n))


if __name__ == "__main__":
    demo1()
