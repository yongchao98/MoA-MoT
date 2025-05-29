import numpy as np

# Coefficients of the polynomial 27p^3 + 10p^2 - 81 = 0
coefficients = [27, 10, 0, -81]

# Find the roots of the polynomial
roots = np.roots(coefficients)

# Filter real roots and round them to 4 decimal places
real_roots = [round(root.real, 4) for root in roots if np.isreal(root)]

# Print the real roots as a comma-separated string
print(",".join(map(str, real_roots)))