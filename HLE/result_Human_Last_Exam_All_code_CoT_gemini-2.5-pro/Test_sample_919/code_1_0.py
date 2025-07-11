import sympy

# Define symbolic variables
mu_0, mu, K_0, a, y, d = sympy.symbols('mu_0 mu K_0 a y d', real=True, positive=True)
i_x = sympy.Symbol('i_x') # Unit vector

# Denominator term from the derivation
denominator = (sympy.cosh(a*d) + (mu_0/mu)*sympy.sinh(a*d))**2

# Numerator term from the derivation
numerator = (mu_0/2) * K_0**2 * sympy.sin(a*y)**2

# The force per unit area vector
force_per_area = (numerator / denominator) * i_x

# To make it look like the options, we can construct the string
# Note: Python's print does not render LaTeX, so we format it as a readable string.
force_expression = f"f/area = (mu_0/2) * (K_0**2 * sin(a*y)**2) / (cosh(a*d) + (mu_0/mu)*sinh(a*d))**2 * i_x"

print("The derived formula for the force per unit area is:")
print(force_expression)

# Let's match with the options.
# Our derived formula has a positive sign, a sin^2 term, and the denominator matches option C.
# The formula is: mu_0/2 * K_0^2 * sin^2(ay) / [cosh(ad) + (mu_0/mu)*sinh(ad)]^2 in the x-direction.
# This corresponds to option C.

print("\nComparing this with the given choices:")
print(f"A. f/area = (K_0^2 * sin(a*y) * cos(a*y)) / [cosh(a*d) + (mu_0/mu)*sinh(a*d)] * i_x")
print(f"B. f/area = (mu_0/2) * (K_0^2 * cos^2(a*y)) / [cosh(a*d) + (mu_0/mu)*sinh(a*d)]^2 * i_x")
print(f"C. f/area = (mu_0/2) * (K_0^2 * sin^2(a*y)) / [cosh(a*d) + (mu_0/mu)*sinh(a*d)]^2 * i_x")
print(f"D. f/area = -(mu_0/2) * (K_0^2 * sin^2(a*y)) / [cosh(a*d) + (mu_0/mu)*sinh(a*d)]^2 * i_x")
print(f"E. f/area = (mu_0/2) * (K_0^2 * sin^2(a*y)) / [cosh(a*d) - (mu_0/mu)*sinh(a*d)]^2 * i_x")
print("\nThe correct option is C.")
