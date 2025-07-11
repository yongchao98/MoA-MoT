import math

# The problem boils down to finding the value of a constant 'k' which represents
# the maximum distance from point A to the surface of the playdough.
# The formula for 'k' is derived from setting the volume of the optimal shape equal to 1.
# Volume V = (4 * pi * k^3) / 15
# Given V = 1 m^3, we solve for k:
# 1 = (4 * pi * k^3) / 15
# 15 = 4 * pi * k^3
# k^3 = 15 / (4 * pi)
# k = (15 / (4 * pi))^(1/3)

# Define the constants from the equation
numerator = 15
denominator_factor = 4

# Calculate the value of k
# This represents the furthest distance in meters.
k_cubed = numerator / (denominator_factor * math.pi)
k = k_cubed**(1/3)

# Print the equation and the final answer
print(f"The equation for the furthest distance (k) is ( {numerator} / ({denominator_factor} * {math.pi}) )^(1/3)")
print(f"The furthest point on the surface of the playdough is {k} meters from point A.")