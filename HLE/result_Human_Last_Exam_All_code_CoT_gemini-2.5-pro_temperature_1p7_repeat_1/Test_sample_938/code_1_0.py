import numpy as np
import math

# Define the coefficients of the characteristic polynomial x^4 - 88x - 57 = 0
# The coefficients are for x^4, x^3, x^2, x^1, x^0
coeffs = [1, 0, 0, -88, -57]

# Find the roots of the polynomial
roots = np.roots(coeffs)

# Calculate the magnitudes of the roots
magnitudes = np.abs(roots)

# Sort the magnitudes in descending order
sorted_magnitudes = sorted(magnitudes, reverse=True)

# The dominant root R has the largest magnitude
R = sorted_magnitudes[0]

# The root with the second largest magnitude is |r_2|
r2_magnitude = sorted_magnitudes[1]

# The problem asks for lim_{n->inf} ln(s_n)/n, which we found to be ln(R * |r_2|)
limit_val = math.log(R * r2_magnitude)

# We need to compute the integer part of 10^4 * limit_val
result = 10000 * limit_val
integer_part = int(result)

# For clarity in the output, let's print the roots and intermediate values
# Find the actual roots corresponding to the magnitudes
R_root = None
r2_root = None
for r in roots:
    if abs(abs(r) - R) < 1e-9:
        R_root = r
    if abs(abs(r) - r2_magnitude) < 1e-9:
        r2_root = r

print(f"The characteristic equation is x^4 - 88x - 57 = 0")
print(f"The roots of the equation are: {roots}")
print(f"The root with the largest magnitude is R = {R_root.real if np.isreal(R_root) else R_root}")
print(f"The root with the second largest magnitude is r_2, with |r_2| = {r2_magnitude}")
print(f"The dominant root magnitude R is: {R}")
print(f"The subdominant root magnitude |r_2| is: {r2_magnitude}")
print(f"The limit L = lim(ln(s_n)/n) is ln(R*|r_2|) = {limit_val}")
print(f"The value to compute is 10^4 * L = {result}")
print(f"The integer part of 10^4 * L is: {integer_part}")