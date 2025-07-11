import math

# This script calculates the maximum distance from point A to the surface of the playdough.
# The shape that maximizes the gravitational field is described by r(θ) = K * sqrt(cos(θ)).
# The constant K is determined by the volume constraint V = 1 m^3.
# The volume of this shape is V = (4 * pi * K^3) / 15.
# Setting V=1 gives K = (15 / (4 * pi))^(1/3).
# The furthest distance from A (the origin) is r_max = r(0) = K.

# Define the constants in the equation K = (numerator / (4 * pi))^(1/3)
numerator = 15
denominator_factor = 4

# Calculate K
pi = math.pi
k = (numerator / (denominator_factor * pi))**(1/3)

# Print the equation and the final answer
print("The equation for the furthest distance (K) is (numerator / (denominator_factor * pi))^(1/3)")
print(f"Substituting the values: K = ({numerator} / ({denominator_factor} * {pi}))^(1/3)")
print(f"The furthest distance is: {k} meters")