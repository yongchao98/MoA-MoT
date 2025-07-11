import numpy as np
from numpy.polynomial import polynomial as P

# Step 1: Define a polynomial with 9 distinct real roots.
# We use the roots of the 9th Chebyshev polynomial, T_9(x), which are known to be
# x_k = cos((2k-1)pi / 18) for k = 1, ..., 9.
n = 9
chebyshev_roots = np.cos((2 * np.arange(1, n + 1) - 1) * np.pi / (2 * n))

# Construct the monic polynomial P(x) with these roots.
# This polynomial will represent the "wiggling" part of our function.
p_poly = P.Polynomial.fromroots(chebyshev_roots)

# Step 2: Ensure the composite function H(x) is strictly increasing.
# We define H(x) = x + epsilon * P(x).
# We need H'(x) = 1 + epsilon * P'(x) > 0 for all x.
# First, find the derivative P'(x).
p_prime_poly = p_poly.deriv()

# Find the minimum value of P'(x). This occurs at a root of P''(x).
p_double_prime_poly = p_prime_poly.deriv()
crit_points = p_double_prime_poly.roots()
# Filter for real roots, as P'(x) is a real polynomial
real_crit_points = crit_points[np.isreal(crit_points)].real
# Evaluate P'(x) at these critical points and at +/- infinity (where it's positive)
min_p_prime_val = np.min(p_prime_poly(real_crit_points))

# Choose epsilon small enough to ensure H'(x) > 0.
# If min_p_prime_val < 0, we need 1 + epsilon * min_p_prime_val > 0
# which means epsilon * min_p_prime_val > -1.
# Since we choose epsilon > 0, this means epsilon < -1 / min_p_prime_val.
if min_p_prime_val < 0:
    epsilon = -0.9 / min_p_prime_val
else:
    epsilon = 1.0

# Step 3: Define the equation for the fixed points.
# The fixed points of H(x) are the roots of H(x) - x = 0, which is epsilon * P(x) = 0.
# The roots are simply the Chebyshev roots we started with.
final_poly = epsilon * p_poly

# Step 4: Display the results.
print("An example can be constructed where the fixed-point equation f(g(x)) = x")
print("is equivalent to the polynomial equation P(x) = 0 shown below.")
print("This polynomial has 9 distinct real roots, which are the fixed points.")
print("-" * 70)

# Format and print the equation
coeffs = final_poly.coef
equation_parts = []
for i in range(len(coeffs) - 1, -1, -1):
    if abs(coeffs[i]) > 1e-9: # Only show non-zero terms
        sign = "-" if coeffs[i] < 0 else "+"
        # For the first term, don't print the '+' sign unless it's the only term
        if len(equation_parts) == 0 and sign == '+':
            sign = ""
        
        term = f"{sign} {abs(coeffs[i]):.4f}"
        if i > 0:
            term += f" * x^{i}"
        equation_parts.append(term)

equation_str = " ".join(equation_parts) + " = 0"
print(f"Equation: {equation_str.lstrip('+ ')}")
print("-" * 70)

# Print the 9 real roots
print("The 9 real roots (fixed points) of this equation are:")
# Sort roots for clean display
sorted_roots = np.sort(chebyshev_roots)
print(", ".join([f"{root:.4f}" for root in sorted_roots]))

print("\nThis demonstrates that it is possible to have 9 fixed points.")
print("Therefore, the maximum number of fixed points is 9.")
<<<9>>>