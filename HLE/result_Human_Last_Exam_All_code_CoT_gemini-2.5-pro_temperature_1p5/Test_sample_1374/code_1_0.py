import math

# The problem requires finding the furthest point 'C' on a specially shaped
# 1 cubic meter volume designed to maximize gravity at a point.

# The volume V of this optimal shape is related to its maximum extent C by the formula:
# V = (4 * pi * C^3) / 15

# We are given V = 1 m^3 and need to solve for C.
# 1 = (4 * pi * C^3) / 15
# Rearranging for C gives:
# C = (15 / (4 * pi))^(1/3)

# Define the numbers in the final equation
volume = 1.0
numerator = 15.0
denominator_constant = 4.0
pi_value = math.pi
exponent = 1/3

# Calculate C
c_cubed = (numerator * volume) / (denominator_constant * pi_value)
max_distance = c_cubed ** exponent

# Print the final equation with its components as requested
print("The final equation for the distance is C = (numerator / (denominator_constant * pi)) ^ exponent")
print("Substituting the values gives:")
print(f"C = ({numerator} / ({denominator_constant} * {pi_value})) ^ ({exponent})")
print(f"\nThe calculation is ({c_cubed}) ^ ({exponent})")
print(f"\nThe furthest point on the surface of the playdough is {max_distance:.9f} meters from point A.")
