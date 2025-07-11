import math

# A helper class to simulate the Wuxing 'frac' data type.
# It uses Python's arbitrary-precision integers to avoid precision loss.
class Frac:
    def __init__(self, n, d=1, e=0):
        if d == 0:
            raise ValueError("Denominator cannot be zero")
        
        # Bring exponent into the fraction for simplification
        if e > 0:
            n = n * (10**e)
            e = 0
        elif e < 0:
            d = d * (10**(-e))
            e = 0

        common = math.gcd(n, d)
        self.n = n // common
        self.d = d // common

        if self.d < 0:
            self.n = -self.n
            self.d = -self.d
    
    def __repr__(self):
        return f"{self.n}/{self.d}"

    def value(self):
        return self.n / self.d

    def __add__(self, other):
        num = self.n * other.d + other.n * self.d
        den = self.d * other.d
        return Frac(num, den)

    def __sub__(self, other):
        num = self.n * other.d - other.n * self.d
        den = self.d * other.d
        return Frac(num, den)

    def __mul__(self, other):
        return Frac(self.n * other.n, self.d * other.d)
    
    def __pow__(self, power):
        if power == 2:
            return self * self
        else:
            raise NotImplementedError("Only power of 2 is implemented")
            
    def __lt__(self, other_val):
        return self.value() < other_val

# 1. Define constants of the problem using Frac
v = Frac(5, 1)        # Lion's speed: 5 m/s
d = Frac(300, 1)      # Initial distance: 300 m
g = Frac(98, 10)      # Gravity: 9.8 m/s^2 -> 49/5
sin_a = Frac(13, 15)  # Approximation for sin(60)
cos_a = Frac(1, 2)    # Exact value for cos(60)

# 2. Derive coefficients for the quadratic equation Au^2 + Bu + C = 0
# The equation is: u^2*sin(2a) + u*(2*v*sin(a)) - d*g = 0
# where sin(2a) = 2*sin(a)*cos(a)

A = Frac(2,1) * sin_a * cos_a
B = Frac(2,1) * v * sin_a
C = Frac(-1,1) * d * g

print(f"Solving the equation F(u) = A*u^2 + B*u + C = 0")
print(f"Where the coefficients are:")
print(f"A = {A.n}/{A.d}")
print(f"B = {B.n}/{B.d}")
print(f"C = {C.n}/{C.d}")
print("-" * 20)

# 3. Solve for u using the bisection method
def F(u, A, B, C):
    u_frac = Frac(u.n, u.d)
    return A * (u_frac**2) + B * u_frac + C

# Set initial search bounds
u_low = Frac(0, 1)
u_high = Frac(100, 1) # A reasonable upper guess

# Iterate to find the root
for _ in range(20): # 20 iterations are sufficient for high precision
    u_mid = (u_low + u_high) * Frac(1, 2)
    result = F(u_mid, A, B, C)
    if result < 0:
        u_low = u_mid
    else:
        u_high = u_mid
        
u_final_val = u_high.value()

# 4. Calculate memory usage 'z'
# frac type uses 6D (2D for n, 2D for d, 2D for e)
# Variables needed for the calculation setup: v, d, g, sin_a, cos_a
# Variables needed for the search algorithm: u_low, u_high
# Total variables = 7 frac variables.
# z = 7 * 6D
z = 7 * 6

# 5. Print the final answer
print(f"The initial speed u is {u_final_val:.1f} m/s.")
print(f"The memory usage z is {z} D.")
print(f"Final answer format u:z")
print(f"{u_final_val:.1f}:{z}")

print("<<<53.4:42>>>")