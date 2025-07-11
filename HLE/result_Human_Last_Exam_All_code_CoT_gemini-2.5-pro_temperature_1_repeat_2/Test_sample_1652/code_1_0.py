import math

# This class simulates the Wuxing architecture's 'frac' data type.
class Frac:
    """
    A custom fraction class to simulate Wuxing's decimal-based computing.
    Value is calculated as (n/d) * 10^e.
    n: signed char (2D, represents -99 to 99)
    d: unsigned char (2D, represents 1 to 99)
    e: signed char (2D, represents -99 to 99)
    """
    def __init__(self, n=0, d=1, e=0):
        if d == 0:
            # In a real scenario, this would be a hardware fault.
            raise ValueError("Denominator cannot be zero")
        self.n = int(n)
        self.d = int(d)
        self.e = int(e)
        self._simplify()

    def _simplify(self):
        """Simulates automatic hardware simplification and truncation."""
        if self.n == 0:
            self.d = 1
            self.e = 0
            return

        # Reduce fraction by greatest common divisor
        common = math.gcd(abs(self.n), self.d)
        self.n //= common
        self.d //= common

        # Truncate numerator to fit in a 2D signed char (-99 to 99)
        while abs(self.n) > 99:
            self.n //= 10
            self.e += 1
        
        # Truncate denominator to fit in a 2D unsigned char (1 to 99)
        while self.d > 99:
            # Round to nearest before truncation
            self.d = (self.d + 5) // 10
            self.e -= 1
        
        if self.d == 0: self.d = 1

    def to_float(self):
        """Converts Frac to a standard float for display purposes."""
        return (self.n / self.d) * (10 ** self.e)

    def __add__(self, other):
        # Align exponents to perform addition with integer math
        common_e = min(self.e, other.e)
        num1 = self.n * (10**(self.e - common_e))
        num2 = other.n * (10**(other.e - common_e))
        new_num = num1 * other.d + num2 * self.d
        new_den = self.d * other.d
        return Frac(new_num, new_den, common_e)

    def __sub__(self, other):
        return self + (other * Frac(-1))

    def __mul__(self, other):
        return Frac(self.n * other.n, self.d * other.d, self.e + other.e)

    def __pow__(self, power):
        if power == 2: return self * self
        raise NotImplementedError("Power only implemented for 2")
    
    def sign(self):
        if self.n > 0: return 1
        if self.n < 0: return -1
        return 0

# --- Main Program ---

# 1. Define constants and variables as Frac objects
# Physical constants
g = Frac(98, 1, -1)      # 9.8
v = Frac(5, 1, 0)        # 5
dist = Frac(3, 1, 2)     # 300

# For angle a = 60 degrees, sin(a) = sin(120) ≈ 0.866
# sin(120) is for sin(2a)
sin60 = Frac(866, 1, -3) # 0.866

# 2. Calculate coefficients for the equation A*u^2 + B*u + C = 0
# A = sin(2a) = sin(120) = sin(60)
A = sin60
# B = 2 * v * sin(a)
B = Frac(2) * v * sin60
# C = -d * g
C = Frac(-1) * dist * g

# As required, output the numbers in the final equation
# We use .to_float() for a readable representation
print(f"The equation to solve is: {A.to_float():.3f} * u^2 + {B.to_float():.3f} * u + {C.to_float():.0f} = 0")

# 3. Solve for u using the Bisection Method (since sqrt is not available)
low = Frac(1, 1, 1)   # Start search at 10 m/s
high = Frac(1, 1, 2)  # End search at 100 m/s

# The function f(u) = Au^2 + Bu + C
def f(u, A, B, C):
    return A * (u**2) + B * u + C

# Iterate 20 times for good precision
for _ in range(20):
    mid = (low + high) * Frac(5, 1, -1) # mid = (low + high) / 2
    if f(mid, A, B, C).sign() == f(low, A, B, C).sign():
        low = mid
    else:
        high = mid

# The result u is the midpoint of the final interval
u_result = mid
u_val = u_result.to_float()

# 4. Calculate memory usage 'z'
# 9 variables (g, v, dist, sin60, B, C, low, high, mid) of type frac (6D each)
z_val = 9 * 6

# 5. Print the final answer in the specified format "u:z"
print(f"\nFinal Answer:")
print(f"{u_val:.1f}:{z_val}")

# <<<53.5:54>>>