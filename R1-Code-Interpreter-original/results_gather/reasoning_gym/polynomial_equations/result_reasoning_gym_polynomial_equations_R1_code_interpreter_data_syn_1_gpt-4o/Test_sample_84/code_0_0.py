import numpy as np

# Coefficients of the polynomial -95z^3 + 62z^2 - 111 = 0
coefficients = [-95, 62, 0, -111]

# Find the roots of the polynomial
roots = np.roots(coefficients)

# Filter real roots and round them to 4 decimal places
real_roots = [root.real for root in roots if np.isclose(root.imag, 0)]
real_roots_rounded = [round(root, 4) for root in real_roots]

# Print the real roots as comma-separated values
print(','.join(map(str, real_roots_rounded)))