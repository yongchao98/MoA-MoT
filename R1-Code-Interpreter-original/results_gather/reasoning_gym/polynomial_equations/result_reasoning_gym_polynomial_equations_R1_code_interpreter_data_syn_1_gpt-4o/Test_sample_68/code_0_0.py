import numpy as np

# Coefficients of the polynomial 85*c^3 - 27*c^2 + 0*c + 40 = 0
coefficients = [85, -27, 0, 40]

# Find the roots of the polynomial
roots = np.roots(coefficients)

# Filter real roots and round them to 4 decimal places
real_roots = [root.real for root in roots if np.isclose(root.imag, 0)]
real_roots_rounded = [round(root, 4) for root in real_roots]

# Print the real roots
print(real_roots_rounded)