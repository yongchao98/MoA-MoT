import math
from fractions import Fraction

def sqrt_fraction(n, iterations=10):
    """
    Calculates the square root of a Fraction using Newton's method.
    This is necessary because the target architecture does not have a built-in sqrt function.
    """
    # An initial guess is required. We use floating point math just for the initial guess,
    # as the subsequent calculations are all done with Fractions.
    x = Fraction(math.sqrt(float(n)))
    for _ in range(iterations):
        if x == 0:
            return Fraction(0)
        # Newton's method iteration: x_new = (x_old + n / x_old) / 2
        x = (x + n / x) / 2
    return x

# --- Main Program ---

# 1. Define problem parameters according to the scenario.
# We use Fractions to adhere to the architecture's no-floating-point rule.
# The rock mass is included as it's part of the problem description.
param_mass = Fraction(1, 2)      # 0.5 kg
param_dist = 300                 # 300 m
param_v = Fraction(5)            # 5 m/s
param_angle = 60                 # 60 degrees
param_g = Fraction(98, 10)       # 9.8 m/s^2

# 2. The physics of the problem leads to a quadratic equation for the initial speed u:
#    (sin(2a)) * u^2 + (2*v*sin(a)) * u - (dist*g) = 0
#    Let's define the coefficients A, B, and C for A*u^2 + B*u + C = 0.

# Calculate trigonometric values using our fraction-based sqrt
# a = 60 deg, 2a = 120 deg
# sin(a) = sin(60) = sqrt(3)/2
# sin(2a) = sin(120) = sqrt(3)/2
sqrt_3 = sqrt_fraction(Fraction(3))
sin_a = sqrt_3 / 2
sin_2a = sin_a

A = sin_2a
B = 2 * param_v * sin_a
C = -param_dist * param_g

# 3. Solve the quadratic equation for u using the formula:
#    u = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
# We take the positive root because speed must be a positive value.
discriminant = B**2 - 4 * A * C
sqrt_of_discriminant = sqrt_fraction(discriminant)

u_frac = (-B + sqrt_of_discriminant) / (2 * A)
u_float = float(u_frac)

# 4. As requested, output each number in the final equation.
print("--- Solving u = (-B + sqrt(B^2 - 4*A*C)) / (2*A) ---")
print(f"B = 2 * v * sin(a) = {float(B):.4f}")
print(f"A = sin(2a) = {float(A):.4f}")
print(f"C = -dist * g = {float(C):.4f}")
print(f"B^2 = {float(B**2):.4f}")
print(f"-4*A*C = {float(-4*A*C):.4f}")
print(f"Discriminant (B^2 - 4AC) = {float(discriminant):.4f}")
print(f"sqrt(Discriminant) = {float(sqrt_of_discriminant):.4f}")
print(f"Numerator (-B + sqrt(Disc.)) = {float(-B + sqrt_of_discriminant):.4f}")
print(f"Denominator (2*A) = {float(2*A):.4f}")
print(f"Final speed u = {u_float:.4f}")
print("-----------------------------------------------------")

# 5. Calculate total memory usage 'z' in decimal positions (D).
# sizeof(int) = 5D; sizeof(frac) = 6D (3 * 2D for char fields)
mem_mass = 6
mem_dist = 5
mem_v = 6
mem_angle = 5
mem_g = 6
mem_u_result = 6
z = mem_mass + mem_dist + mem_v + mem_angle + mem_g + mem_u_result

# 6. Format the final answer as u:z, with u rounded to one decimal place.
u_rounded = round(u_float, 1)
print(f"The final answer is u:{u_rounded} and z:{z}")
<<<53.5:34>>>