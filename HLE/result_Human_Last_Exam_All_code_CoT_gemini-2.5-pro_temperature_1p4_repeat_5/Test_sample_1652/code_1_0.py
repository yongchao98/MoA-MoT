import math

# This script simulates calculations on the WUXING architecture to find the
# required initial speed 'u' for the rock.

def gcd(a, b):
    """Helper function to find the greatest common divisor."""
    return math.gcd(a, b)

class Frac:
    """
    A class to simulate the WUXING 'frac' type.
    Represents a number as (numerator/denominator) * 10^exponent.
    """
    def __init__(self, n, d=1, e=0):
        if d == 0:
            raise ZeroDivisionError
        self.n = n
        self.d = d
        self.e = e
        self._simplify()

    def _simplify(self):
        """Simplifies the fraction and normalizes the exponent."""
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d
        if self.n == 0:
            self.d = 1
            self.e = 0
            return
        common = gcd(abs(self.n), self.d)
        self.n //= common
        self.d //= common
        while self.n != 0 and self.n % 10 == 0:
            self.n //= 10
            self.e += 1
        while self.d != 0 and self.d % 10 == 0:
            self.d //= 10
            self.e -= 1

    def to_float(self):
        """Converts the fraction to a floating-point number for final display."""
        return (self.n / self.d) * (10 ** self.e)

    def __repr__(self):
        """String representation for debugging."""
        return f"{self.n}/{self.d}e{self.e}"

    def __add__(self, other):
        e_new = max(self.e, other.e)
        n1_adj = self.n * (10 ** (e_new - self.e))
        n2_adj = other.n * (10 ** (e_new - other.e))
        new_n = n1_adj * other.d + n2_adj * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d, e_new)

    def __sub__(self, other):
        return self + Frac(-other.n, other.d, other.e)

    def __mul__(self, other):
        return Frac(self.n * other.n, self.d * other.d, self.e + other.e)

    def __truediv__(self, other):
        return Frac(self.n * other.d, self.d * other.n, self.e - other.e)

# 1. Define problem parameters as Frac objects
dist = Frac(300)       # Initial distance = 300 m
v = Frac(5)            # Lion speed v = 5 m/s
g = Frac(98, 1, -1)    # Gravity g = 9.8 m/s^2
sin_a = Frac(433, 500) # sin(60) approximated as 0.866
two = Frac(2)

# 2. Formulate the equation: a*u^2 + b*u + c = 0
# Simplified to: 1*u^2 + (2*v)*u - (dist*g/sin_a) = 0
A_coeff = Frac(1)
B_coeff = two * v
C_coeff = (dist * g / sin_a) * Frac(-1) # The term is negative

print("The problem is solved using the equation: a*u^2 + b*u + c = 0")
print(f"Coefficient a: {A_coeff.to_float()}")
print(f"Coefficient b: {B_coeff.to_float()}")
print(f"Coefficient c: {C_coeff.to_float()}")
print("---")

# 3. Solve for 'u' using Newton-Raphson iteration
# f(u) = A*u^2 + B*u + C
# f'(u) = 2*A*u + B
# u_new = u - f(u)/f'(u)
u = Frac(50)  # Initial guess u = 50 m/s

for _ in range(10): # 10 iterations is sufficient for convergence
    f_u = (A_coeff * u * u) + (B_coeff * u) + C_coeff
    f_prime_u = (two * A_coeff * u) + B_coeff
    if f_prime_u.n == 0:
        break
    u = u - (f_u / f_prime_u)

# 4. Final result for 'u'
u_result = u.to_float()

# 5. Calculate memory usage 'z'
# 5 main variables (u, v, dist, g, sin_a), each a 'frac' of size 6D.
z = 5 * 6

# 6. Print the answer in the specified format u:z
print(f"Final Answer (u:z):")
print(f"{u_result:.1f}:{z}")