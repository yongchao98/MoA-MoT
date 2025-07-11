import numpy as np

# The problem is reduced to finding the real root of a cubic equation for the position x0.
# The equation is: x^3 - 6x^2 + 18x - 18 = 0.

# Define the coefficients of the polynomial.
a = 1
b = -6
c = 18
d = -18
coeffs = [a, b, c, d]

# Print the numbers (coefficients) of the final equation being solved.
print(f"The final equation for x0 is a cubic polynomial: {a}x^3 + ({b})x^2 + {c}x + ({d}) = 0")

# Find the roots of the polynomial using numpy.
roots = np.roots(coeffs)

# The solution must be a real number, as it represents a physical position.
# We iterate through the roots to find the real one.
real_root = None
for root in roots:
    if np.isreal(root):
        real_root = np.real(root)
        break

# Print the final result.
print("\nThe position x0 where the particle reaches y = -3 is:")
print(real_root)
