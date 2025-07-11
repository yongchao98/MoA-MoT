import math

# This script simulates calculations on the Wuxing architecture to find the rock's initial speed.

def gcd(a, b):
    """Helper function to compute the greatest common divisor."""
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a

class Frac:
    """
    A class to simulate the 'frac' data type of the Wuxing architecture.
    It represents a number as (n/d) * 10^e.
    """
    def __init__(self, n, d=1, e=0):
        if d == 0:
            raise ZeroDivisionError
        # Use large integers to maintain precision, simulating the Wuxing compiler's features.
        self.n = int(n)
        self.d = int(d)
        self.e = int(e)
        self._simplify()

    def _simplify(self):
        if self.n == 0:
            self.d = 1
            return
        common = gcd(self.n, self.d)
        self.n //= common
        self.d //= common
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d

    def to_float(self):
        """Converts the fraction to a standard float for display."""
        return (self.n / self.d) * (10 ** self.e)

    def __add__(self, other):
        # Align exponents: a/b * 10^e1 + c/d * 10^e2
        if self.e > other.e:
            num2 = other.n * (10**(self.e - other.e))
            den2 = other.d
            new_e = self.e
            num1, den1 = self.n, self.d
        else:
            num1 = self.n * (10**(other.e - self.e))
            den1 = self.d
            new_e = other.e
            num2, den2 = other.n, other.d
        
        num = num1 * den2 + num2 * den1
        den = den1 * den2
        return Frac(num, den, new_e)

    def __neg__(self):
        return Frac(-self.n, self.d, self.e)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return Frac(self.n * other.n, self.d * other.d, self.e + other.e)

    def __truediv__(self, other):
        return Frac(self.n * other.d, self.d * other.n, self.e - other.e)

def sqrt_frac(s_frac, iterations=20):
    """
    Calculates the square root of a Frac using the Babylonian method.
    This is required as sqrt() is not a built-in function.
    """
    if s_frac.n < 0:
        raise ValueError("Cannot take sqrt of a negative number.")
    if s_frac.n == 0:
        return Frac(0)

    # Adjust exponent to be even for sqrt(10^e)
    if s_frac.e % 2 != 0:
        s_frac.n *= 10
        s_frac.e -= 1
    
    new_e = s_frac.e // 2
    val_to_sqrt = Frac(s_frac.n, s_frac.d)

    # Initial guess to speed up convergence
    guess_f = math.sqrt(val_to_sqrt.to_float())
    x = Frac(int(guess_f * 10**8), 10**8)
    
    two = Frac(2)
    for _ in range(iterations):
        x = (x + val_to_sqrt / x) / two
        
    x.e += new_e
    return x

# --- Main Program ---

# 1. Define problem variables based on the Wuxing data types
v = Frac(5)          # lion's speed (m/s)
g = Frac(98, 10)     # gravity (m/s^2)
distance = Frac(300) # initial distance (m)
angle_a = 60         # launch angle (degrees)

# The derived quadratic equation for the initial speed u is:
# u^2 + (2*v)*u - (distance * g / sin(angle_a)) = 0
# This is in the form A*u^2 + B*u + C = 0

# 2. Calculate coefficients A, B, and C using fractional arithmetic
# We need sin(60) = sqrt(3)/2
three = Frac(3)
two = Frac(2)
sqrt3 = sqrt_frac(three)
sin60 = sqrt3 / two

A_coeff = Frac(1)
B_coeff = two * v
C_coeff = -(distance * g / sin60)

# 3. Solve for u using the quadratic formula: u = (-B + sqrt(B^2 - 4AC)) / 2A
four = Frac(4)
discriminant = (B_coeff * B_coeff) - (four * A_coeff * C_coeff)
sqrt_discriminant = sqrt_frac(discriminant)
numerator = -B_coeff + sqrt_discriminant
denominator = two * A_coeff
u_final = numerator / denominator

# 4. Print the breakdown of the final equation as requested
print("Solving u = (-B + sqrt(B^2 - 4*A*C)) / (2*A)")
print(f"B = {B_coeff.to_float()}")
print(f"B^2 = {(B_coeff * B_coeff).to_float()}")
print(f"4*A*C = {(four * A_coeff * C_coeff).to_float()}")
print(f"B^2 - 4*A*C = {discriminant.to_float()}")
print(f"sqrt(B^2 - 4*A*C) = {sqrt_discriminant.to_float()}")
print(f"Numerator (-B + sqrt(...)) = {numerator.to_float()}")
print(f"Denominator (2*A) = {denominator.to_float()}")
print(f"Final u = {u_final.to_float()}")
print("-" * 20)

# 5. Calculate memory usage 'z'
# Variables: v, g, distance, m (mass, unused but given), u -> 5 frac
#            angle_a -> 1 char
# size_frac = 3 * size_char = 3 * 2D = 6D
# size_char = 2D
# z = 5 * 6D + 1 * 2D = 32D
z = 32

# 6. Format and print the final answer
u_rounded = round(u_final.to_float(), 1)
print(f"{u_rounded}:{z}")