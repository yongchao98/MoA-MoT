import numpy as np

# Coefficients of the polynomial 97*m^3 - 66*m^2 + 25*m - 34
coefficients = [97, -66, 25, -34]

# Find the roots of the polynomial
roots = np.roots(coefficients)

# Filter out the real roots and round them to 4 decimal places
real_roots = [root.real for root in roots if np.isclose(root.imag, 0)]
real_roots_rounded = [round(root, 4) for root in real_roots]

# Print the real roots as a comma-separated string
print(",".join(map(str, real_roots_rounded)))