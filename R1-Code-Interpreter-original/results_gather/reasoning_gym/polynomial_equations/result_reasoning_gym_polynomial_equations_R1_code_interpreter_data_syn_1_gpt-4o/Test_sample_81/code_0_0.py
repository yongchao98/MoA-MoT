import numpy as np

# Coefficients of the polynomial -148a^3 + 56a^2 - 9 = 0
coefficients = [-148, 56, 0, -9]

# Find the roots of the polynomial
roots = np.roots(coefficients)

# Filter real roots and round them to 4 decimal places
real_roots = [round(root.real, 4) for root in roots if np.isreal(root)]

# Print the real roots as comma-separated values
print(",".join(map(str, real_roots)))