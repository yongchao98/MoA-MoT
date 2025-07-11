import math

# This script calculates the furthest distance from a point A to the surface of a 
# 1 cubic meter volume arranged to create the strongest possible gravitational field at A.

# As derived from the physics principles, the shape of the object is described by
# r(θ) = k * sqrt(cos(θ)).
# The volume V=1 is related to k by the equation: V = (4 * pi * k^3) / 15.
# The furthest point on the surface is at a distance r_max = k.
# We solve for k: k = (15 / (4 * pi))^(1/3).

# Define the numbers in the final equation
numerator = 15
denominator_part = 4
power = 1/3

# We will use the math.pi constant for pi.
print("The final equation for the furthest distance (d) is:")
print(f"d = ({numerator} / ({denominator_part} * π)) ^ ({power})")
print("-" * 30)

# Calculate the numerical value
distance = (numerator / (denominator_part * math.pi)) ** power

print("The furthest distance from A in meters is:")
print(distance)
