import math

class Frac:
    """
    A class to simulate the Wuxing 'frac' type and its arithmetic.
    It includes automatic simplification of powers of 10.
    """
    def __init__(self, n, d=1, e=0):
        if not all(isinstance(i, int) for i in [n, d, e]):
            raise TypeError("Numerator, denominator, and exponent must be integers.")
        if d == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        self.n = n
        self.d = d
        self.e = e
        self.simplify()

    def simplify(self):
        """Simplifies the fraction and adjusts the exponent."""
        if self.n == 0:
            self.d = 1
            self.e = 0
            return

        # Simulate Wuxing simplification: pull factors of 10 from n/d to exponent
        while self.n % 10 == 0 and self.n != 0:
            self.n //= 10
            self.e += 1
        while self.d % 10 == 0 and self.d != 0:
            self.d //= 10
            self.e -= 1

        # Standard simplification by GCD
        common_divisor = math.gcd(self.n, self.d)
        self.n //= common_divisor
        self.d //= common_divisor

        # Ensure denominator is positive
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d

    def __mul__(self, other):
        if isinstance(other, int):
            other = Frac(other)
        return Frac(self.n * other.n, self.d * other.d, self.e + other.e)

    def __truediv__(self, other):
        if isinstance(other, int):
            other = Frac(other)
        return Frac(self.n * other.d, self.d * other.n, self.e - other.e)

    def __sub__(self, other):
        if isinstance(other, int):
            other = Frac(other)
        # Find common exponent
        new_e = max(self.e, other.e)
        n1 = self.n * (10**(self.e - new_e))
        n2 = other.n * (10**(other.e - new_e))
        
        # Perform subtraction with common denominator
        common_d = self.d * other.d
        new_n = n1 * other.d - n2 * self.d
        
        return Frac(int(new_n), common_d, new_e)

    def __rsub__(self, other):
        # Handles cases like `1 - frac`
        if isinstance(other, int):
            other = Frac(other)
        return other - self

    def __repr__(self):
        return f"{self.n}/{self.d}e{self.e}"

    @property
    def value(self):
        return (self.n / self.d) * (10**self.e)

# 1. Define constants for the problem
# G ≈ 6.67e-11 ≈ 20/3e-11
G = Frac(20, 3, -11)
# M = 2 * M_sun ≈ 2 * 1.989e30 kg ≈ 4e30 kg
M = Frac(4, 1, 30)
# r = R + d = 20km + 60km = 80km = 8e4 m
r = Frac(8, 1, 4)
# c ≈ 3e8 m/s
c = Frac(3, 1, 8)
two = Frac(2)
one = Frac(1)

# 2. Calculate x = (2*G*M) / (r*c^2)
print("Calculating f using the approximation f ≈ 1 - x/2, where x = (2*G*M)/(r*c^2)\n")

# Numerator
numer = two * G * M
# Denominator
denom = r * c * c
# x
x = numer / denom

print("Equation steps:")
print(f"  2 * G * M = {repr(numer)}")
print(f"  r * c^2   = {repr(denom)}")
print(f"  x = {repr(numer)} / {repr(denom)} = {repr(x)}")

# 3. Calculate f ≈ 1 - x/2
x_div_2 = x / two
f = one - x_div_2

print(f"  x / 2     = {repr(x)} / 2 = {repr(x_div_2)}")
print(f"  f ≈ 1 - (x / 2) = 1 - {repr(x_div_2)} = {repr(f)}\n")

# 4. Calculate f_val and z_val
f_val = f.value
# Memory usage z:
# struct frac is 6D (2D+2D+2D).
# An efficient C program needs variables for G, M, r, c, plus two temporaries.
# Total = 6 variables * 6D/variable = 36D.
z_val = "36D"

print(f"Final calculated dilation factor f: {f_val:.4f}")
print(f"Estimated memory usage z: {z_val}")
print("---")
print(f"Final Answer f:z = {f_val:.4f}:{z_val}")
