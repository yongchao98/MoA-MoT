import math

# A helper function for fraction simplification
def gcd(a, b):
    return math.gcd(a, b)

class Frac:
    """
    A class to simulate the WUXING frac type and its arithmetic.
    Value is (n/d). Exponent 'e' is handled during initialization.
    """
    def __init__(self, n, d=1, e=0):
        if d == 0:
            raise ZeroDivisionError
        
        # Handle decimal exponent 'e' by modifying n and d
        if e > 0:
            n *= (10**e)
        elif e < 0:
            d *= (10**(-e))
            
        common = gcd(n, d)
        self.n = n // common
        self.d = d // common

        # Convention: keep the denominator positive
        if self.d < 0:
            self.n = -self.n
            self.d = -self.d

    def __repr__(self):
        """String representation for printing."""
        if self.d == 1:
            return f"{self.n}"
        return f"{self.n}/{self.d}"

    def to_float(self):
        """Convert fraction to a standard float for the final result."""
        return self.n / self.d

    def __add__(self, other):
        new_n = self.n * other.d + other.n * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __sub__(self, other):
        new_n = self.n * other.d - other.n * self.d
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __mul__(self, other):
        new_n = self.n * other.n
        new_d = self.d * other.d
        return Frac(new_n, new_d)

    def __truediv__(self, other):
        new_n = self.n * other.d
        new_d = self.d * other.n
        return Frac(new_n, new_d)

# --- Main Program Logic ---

# 1. Define persistent variables for the problem.
v = Frac(5)                   # Lion's speed in m/s
initial_distance = Frac(300)  # Initial distance in meters
g = Frac(98, 10)              # Acceleration due to gravity (9.8 m/s^2)

# Use a rational approximation for sin(60) that is accurate and fits 'char' constraints.
# sin(60) = sqrt(3)/2 ≈ 0.866... The fraction 13/15 ≈ 0.8666... is a good choice.
sin_60_approx = Frac(13, 15)
# sin(120) = sin(60), so we use the same value.

# 2. Set up the quadratic equation: A_eq*u^2 + B_eq*u + C_eq = 0
# derived from u^2*sin(2a) + u*(2*v*sin(a)) - d*g = 0
A_eq = sin_60_approx
B_eq = Frac(2) * v * sin_60_approx
C_eq = initial_distance * g * Frac(-1)

print("The problem reduces to solving the quadratic equation A*u^2 + B*u + C = 0 for u.")
print("The calculated coefficients are:")
print(f"A = {A_eq}")
print(f"B = {B_eq}")
print(f"C = {C_eq}")
print("-" * 20)

# 3. Solve the equation numerically using the Newton-Raphson method.
# f(u) = A_eq*u^2 + B_eq*u + C_eq
# f'(u) = 2*A_eq*u + B_eq
def f(u_frac):
    return A_eq * u_frac * u_frac + B_eq * u_frac + C_eq

def f_prime(u_frac):
    return Frac(2) * A_eq * u_frac + B_eq

# Make an initial guess for u. A rough estimate suggests u is around 50-60.
u = Frac(50) 
# Iterate to find a precise root.
for _ in range(10):
    u = u - f(u) / f_prime(u)

# 4. Finalize results.
# Round the final velocity to one decimal place.
u_final = round(u.to_float(), 1)

# Calculate memory usage 'z' in 'D' (decimal digits).
# There are 4 persistent variables: v, initial_distance, g, sin_60_approx.
# Each 'frac' struct has 3 'char's (n, d, e), and each char is 2D.
# z = 4 vars * 3 components/var * 2D/component.
z = 4 * 3 * 2

print(f"The calculated initial speed u is: {u_final} m/s")
print(f"The memory usage z for stored variables is: {z}D")
print("-" * 20)
print("Final answer in the format u:z")
print(f"{u_final}:{z}")