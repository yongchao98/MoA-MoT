import math

# Helper function for Greatest Common Divisor
def gcd(a, b):
    return math.gcd(a, b)

class Frac:
    """
    A class to simulate the Wuxing 'frac' data type.
    Value is calculated as (n/d) * 10^e.
    n, d, e are constrained by 'char' type in the architecture, but
    this simulation uses Python's arbitrary precision integers for intermediate steps.
    """
    def __init__(self, n, d=1, e=0):
        if not isinstance(n, int) or not isinstance(d, int) or not isinstance(e, int):
            raise TypeError("Numerator, denominator, and exponent must be integers.")
        if d == 0:
            raise ZeroDivisionError("Denominator cannot be zero.")
        
        # Handle conversion from large integers by adjusting the exponent
        if isinstance(n, int) and d == 1 and e == 0:
            val = n
            exp = 0
            while val % 10 == 0 and val != 0:
                val //= 10
                exp += 1
            self.n = val
            self.d = d
            self.e = exp
        else:
            self.n = n
            self.d = d
            self.e = e
        
        self.simplify()

    def simplify(self):
        # Ensure sign is only in the numerator
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d
        
        # Reduce the fraction by the greatest common divisor
        if self.n == 0:
            self.d = 1
        else:
            common = gcd(abs(self.n), self.d)
            self.n //= common
            self.d //= common

    def to_float(self):
        return (self.n / self.d) * (10 ** self.e)

    def __repr__(self):
        return f"({self.n}/{self.d})e{self.e}"

    def __add__(self, other):
        # Align exponents
        max_e = max(self.e, other.e)
        n1 = self.n * (10**(max_e - self.e))
        n2 = other.n * (10**(max_e - other.e))
        
        # Common denominator addition
        new_n = n1 * other.d + n2 * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d, max_e)

    def __sub__(self, other):
        # Implemented as addition with a negative
        neg_other = Frac(-other.n, other.d, other.e)
        return self.__add__(neg_other)

    def __mul__(self, other):
        new_n = self.n * other.n
        new_d = self.d * other.d
        new_e = self.e + other.e
        return Frac(new_n, new_d, new_e)

    def __truediv__(self, other):
        if other.n == 0:
            raise ZeroDivisionError("Fraction division by zero.")
        new_n = self.n * other.d
        new_d = self.d * other.n
        new_e = self.e - other.e
        return Frac(new_n, new_d, new_e)

def frac_sqrt(s, iterations=20):
    """
    Calculates the square root of a Frac using the Babylonian method.
    """
    # Use a floating point representation for a good initial guess
    s_float = s.to_float()
    if s_float < 0:
        raise ValueError("Cannot compute square root of a negative number")
    
    # Start with a reasonable guess
    x = Frac(int(math.sqrt(s_float))) if s_float > 0 else Frac(0)
    
    two = Frac(2)
    for _ in range(iterations):
        if x.n == 0: # Avoid division by zero if guess is 0
            return Frac(0)
        x = (x + s / x) / two
    return x

# Main execution
# 1. Define problem constants according to Wuxing constraints
# angle a = 60 degrees.
# sin(60) = sqrt(3)/2. We approximate sqrt(3) with 26/15, so n,d fit in a char.
# sin(60) approx (26/15)/2 = 13/15.
# sin(2*60) = sin(120) = sin(60), so same approximation.
v_lion = Frac(5)
d_initial = Frac(300)
g = Frac(98, 10)  # 9.8 = 98/10, simplified to 49/5
sin_a = Frac(13, 15)
sin_2a = Frac(13, 15)
TWO = Frac(2)

# 2. Physics: (sin(2a))u^2 + (2*v*sin(a))u - (d*g) = 0
# This is a quadratic equation Au^2 + Bu + C = 0
A = sin_2a
B = TWO * v_lion * sin_a
C = Frac(-1) * d_initial * g

# 3. Solve using the quadratic formula: u = (-B + sqrt(B^2 - 4AC)) / 2A
B_squared = B * B
FOUR_A_C = Frac(4) * A * C
discriminant = B_squared - FOUR_A_C

sqrt_discriminant = frac_sqrt(discriminant)
numerator = (Frac(-1) * B) + sqrt_discriminant
denominator = TWO * A
u_frac = numerator / denominator
u_final = u_frac.to_float()

# 4. Calculate memory usage 'z'
# d_initial (300) fits in a 5D 'int'
# v_lion (5) fits in a 2D 'signed char'
# g is a 'frac' (n: signed char 2D, d: unsigned char 2D, e: signed char 2D) = 6D
# sin_a is a 'frac' = 6D
# sin_2a is a 'frac' = 6D
# z = 5 + 2 + 6 + 6 + 6 = 25 D
z = 25

# 5. Print the results as requested
print("The problem is to find the initial speed 'u' of a projectile to hit a moving target.")
print("The physics of the problem leads to a quadratic equation: A*u^2 + B*u + C = 0")
print("Where A = sin(2a), B = 2*v*sin(a), and C = -d*g\n")

print("--- Calculated Values (Wuxing 'frac' representation) ---")
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}\n")

print("--- Solving the Equation ---")
print("We use the formula: u = (-B + sqrt(B^2 - 4AC)) / 2A")
print(f"B^2 = {B_squared}")
print(f"4AC = {FOUR_A_C}")
print(f"Discriminant (B^2 - 4AC) = {discriminant}")
print(f"Square root of Discriminant = {sqrt_discriminant}")
print(f"Numerator (-B + sqrt(...)) = {numerator}")
print(f"Denominator (2A) = {denominator}")
print(f"u (as fraction) = {u_frac}\n")

print("--- Final Answer ---")
print(f"The initial speed u = {u_final:.1f} m/s")
print(f"The memory usage z = {z} D")

# Format final answer as requested by the user
final_answer = f"{u_final:.1f}:{z}"
print(f"<<<{final_answer}>>>")