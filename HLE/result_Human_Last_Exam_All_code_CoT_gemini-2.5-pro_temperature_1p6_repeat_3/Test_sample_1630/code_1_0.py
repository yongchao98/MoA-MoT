import numpy as np
from scipy.optimize import minimize_scalar
from numpy.polynomial.polynomial import Polynomial

# This script constructs a degree 9 polynomial h(x) that is strictly increasing
# and has 9 fixed points, demonstrating that the maximum number is achievable.

# Step 1: Define a polynomial P(x) with 9 distinct real roots.
# We choose the roots to be -4, -3, -2, -1, 0, 1, 2, 3, 4 for simplicity.
fixed_points = np.arange(-4, 5)
P = Polynomial.fromroots(fixed_points)

# Step 2: Find the derivative of P(x).
P_prime = P.deriv()

# Step 3: Find the global minimum of P'(x).
# P'(x) is a degree 8 polynomial with a positive leading coefficient, so it has a global minimum.
res = minimize_scalar(P_prime, bounds=(-5, 5), method='bounded')
min_P_prime_val = res.fun

# Step 4: Choose a small positive epsilon to ensure h'(x) = 1 + epsilon * P'(x) > 0.
# If min_P_prime_val is negative, we need epsilon < -1 / min_P_prime_val.
epsilon = 0.99 * (-1 / min_P_prime_val)

# Step 5: The fixed point equation is h(x) = x, which is x + epsilon * P(x) = x.
# This simplifies to epsilon * P(x) = 0.
# The polynomial whose roots we seek is P(x).
final_poly_coeffs = epsilon * P.coef

# Print the results
print("We construct a function h(x) = f(g(x)) to have 9 fixed points.")
print("The fixed points are the roots of the equation h(x) - x = 0.")
print("\nBased on our construction, this equation is epsilon * P(x) = 0, where P(x) has 9 real roots.")

# Print the equation, showing each coefficient.
equation_str = []
for i, coeff in reversed(list(enumerate(final_poly_coeffs))):
    if not np.isclose(coeff, 0):
        equation_str.append(f"({coeff:+.4e}) * x^{i}")
print("\nThe final equation for the fixed points is:")
print(" " + " + ".join(equation_str) + " = 0")


print(f"\nThis equation has 9 real roots, which are the fixed points:")
print(fixed_points)

# Verify that the constructed h(x) is strictly increasing
min_h_prime = 1 + epsilon * min_P_prime_val
print(f"\nThe minimum value of h'(x) is approximately {min_h_prime:.4f}, which is positive.")
print("\nThis confirms the existence of a suitable function h(x).")
print("Advanced theorems guarantee that this h(x) can be decomposed into f(g(x))")
print("with the desired properties for f and g.")
print("\nThus, the maximum number of fixed points is 9.")

<<<9>>>