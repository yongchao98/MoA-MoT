import math

# This script calculates the distance to the furthest point on a playdough shape
# designed to create the strongest possible gravitational field at a point A.

# The volume of the playdough is given.
volume = 1.0  # in cubic meters

# The formula for the maximum distance (r_max) is derived from maximizing the
# gravitational field integral, which results in the relationship:
# r_max = (15 * V / (4 * pi))^(1/3)

# We will now calculate this value step-by-step.

# 1. Define the numerator of the fraction
numerator = 15 * volume

# 2. Define the denominator of the fraction
pi_val = math.pi
denominator = 4 * pi_val

# 3. Calculate the ratio
ratio = numerator / denominator

# 4. Calculate the cube root of the ratio to find the final answer
r_max = ratio**(1/3)

# Print the equation and the values used at each step
print("The formula for the furthest distance (r_max) is:")
print(f"r_max = (15 * Volume / (4 * pi))^(1/3)")
print("\nSubstituting the values:")
print(f"r_max = ({numerator} / (4 * {pi_val}))^(1/3)")
print(f"r_max = ({numerator} / {denominator})^(1/3)")
print(f"r_max = ({ratio})^(1/3)")
print("\nThe final result is:")
print(f"r_max = {r_max} meters")
