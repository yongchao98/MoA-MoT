import numpy as np

# Coefficients of the polynomial -29g^3 + 45g^2 - 88g + 99
coefficients = [-29, 45, -88, 99]

# Find the roots of the polynomial
roots = np.roots(coefficients)

# Filter real roots and round them to 4 decimal places
real_roots = [round(root.real, 4) for root in roots if np.isreal(root)]

# Print the real roots as a comma-separated string
print(",".join(map(str, real_roots)))