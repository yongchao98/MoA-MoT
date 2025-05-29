import numpy as np

# Coefficients of the polynomial -46r^3 - 51r^2 + 0*r + 101
coefficients = [-46, -51, 0, 101]

# Find the roots of the polynomial
roots = np.roots(coefficients)

# Filter real roots and round them to 4 decimal places
real_roots = [round(root.real, 4) for root in roots if np.isclose(root.imag, 0)]

# Print the real roots as a comma-separated string
print(",".join(map(str, real_roots)))