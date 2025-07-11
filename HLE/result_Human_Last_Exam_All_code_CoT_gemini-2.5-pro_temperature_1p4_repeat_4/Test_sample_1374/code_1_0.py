import math

# This script calculates the maximum distance from a point A to the surface of a
# 1 cubic meter volume arranged to create the strongest possible gravitational field at A.

# The equation for the distance is derived from the calculus of variations.
# The numbers used in the final equation are:
num = 15
den_1 = 4
den_2_name = "pi"
den_2_val = math.pi
pwr_num = 1
pwr_den = 3
pwr = pwr_num / pwr_den

# Calculate the result
base = num / (den_1 * den_2_val)
distance = base ** pwr

# Print the final equation with all its numerical components
print("The furthest distance 'd' is calculated by the following equation:")
print(f"d = ({num} / ({den_1} * {den_2_name})) ^ ({pwr_num}/{pwr_den})")
print("\nCalculating the result...")
print(f"The furthest point on the surface of the playdough is {distance} meters from point A.")
