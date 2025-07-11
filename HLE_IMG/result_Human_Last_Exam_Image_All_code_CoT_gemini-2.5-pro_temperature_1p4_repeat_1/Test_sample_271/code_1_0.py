# Define the given physical quantities
m1 = 1  # mass at the end of the rod in kg
m2 = 2  # mass of the projectile in kg
R = 3   # length of the rod in m
d = 1   # distance from pivot O to the collision point in m
g = 10  # acceleration due to gravity in m/s^2

# Based on the principles of conservation of angular momentum and energy,
# the following equation for v0^2 can be derived. The code constructs a
# string to display this equation with all the given numbers substituted,
# as requested.
# The general formula is:
# v0^2 = (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2*d)^2

print("The proposed equation for the value that v0 must have is:")
print(f"v0^2 = (4 * {g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}^2 + {m2} * {d}^2)) / ({m2} * {d})^2")