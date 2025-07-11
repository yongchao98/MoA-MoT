import math

# The volume of the playdough is given in cubic meters.
V = 1.0

# Constants from the derived formula r_max = (15 * V / (4 * pi))^(1/3)
numerator_constant = 15.0
denominator_constant = 4.0

# Calculate the maximum distance (r_max)
# This is the distance from point A to the furthest point on the surface.
distance = (numerator_constant * V / (denominator_constant * math.pi))**(1/3)

# The final equation is: distance = (15 * 1 / (4 * pi))^(1/3)
# We print the components and the final result.
print("The calculation for the furthest distance is based on the formula:")
print("distance = (15 * V / (4 * pi))^(1/3)")
print("\nSubstituting the given values:")
print(f"distance = ({numerator_constant} * {V} / ({denominator_constant} * {math.pi}))^(1/3)")
print("\nFinal Result:")
print(f"The furthest point on the surface of the playdough is {distance} meters from point A.")