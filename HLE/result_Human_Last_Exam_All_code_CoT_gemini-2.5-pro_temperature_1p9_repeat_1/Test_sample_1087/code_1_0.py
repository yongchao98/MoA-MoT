import math

# The largest value of r is the diagonal length 'd' of the largest regular
# pentagon that can be inscribed in a unit square.
# This value is calculated as d = sin(72°) / cos(18°).

# Define angles in radians
angle_72_deg_in_rad = math.radians(72)
angle_18_deg_in_rad = math.radians(18)

# Get the values for the numerator and denominator
sin_72 = math.sin(angle_72_deg_in_rad)
cos_18 = math.cos(angle_18_deg_in_rad)

# Calculate the final result for r
r = sin_72 / cos_18

# As requested, here are the numbers from the final equation for r
print("The final equation is of the form: r = sin(72°) / cos(18°)")
print(f"The value for sin(72°) is: {sin_72}")
print(f"The value for cos(18°) is: {cos_18}")
print(f"The final value for r is: {r}")
