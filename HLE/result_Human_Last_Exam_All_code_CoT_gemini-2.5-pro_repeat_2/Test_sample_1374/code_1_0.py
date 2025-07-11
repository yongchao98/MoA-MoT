import math

# The problem is to find the furthest point on a 1 m^3 shape
# that maximizes the gravitational field at a point A.
# The shape is described by r(theta) = k * sqrt(cos(theta)).
# The furthest point is at a distance k from A, where k is determined
# by the volume constraint.

# The equation relating the volume (V=1) to the distance (k) is:
# k^3 = 15 / (4 * pi)
# Therefore, k = (15 / (4 * pi))^(1/3)

# Define the numbers in the final equation
numerator = 15
denominator_coefficient = 4
pi_value = math.pi
power = 1/3

# Calculate the result
distance = (numerator / (denominator_coefficient * pi_value))**power

# As requested, output the numbers in the final equation and the result.
# The equation for the furthest distance (r_max) is derived from the volume constraint.
print(f"The equation for the furthest distance is: ( {numerator} / ({denominator_coefficient} * π) )^(1/3)")
print(f"The calculated distance is: {distance} meters")