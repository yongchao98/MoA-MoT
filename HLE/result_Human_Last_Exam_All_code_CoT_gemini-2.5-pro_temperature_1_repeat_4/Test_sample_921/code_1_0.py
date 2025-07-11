import math
import sys

# Define the input values provided in the problem
angle_in_degrees = 40
distance = 100

# 1. Convert the angle from degrees to radians
# Trigonometric functions in most programming languages require radians
angle_in_radians = math.radians(angle_in_degrees)

# 2. Calculate the tangent of the angle
# This is the ratio of the opposite side (height) to the adjacent side (distance)
tan_of_angle = math.tan(angle_in_radians)

# 3. Calculate the final height
height = distance * tan_of_angle

# As requested, here are the numbers that form the final equation:
# height = distance * tan(radians(angle_in_degrees))
print("--- Calculation Breakdown ---")
print(f"Distance to Palace: {distance} meters")
print(f"Angle to Top: {angle_in_degrees} degrees")
print(f"Resulting Height: {height} meters")
print("\n--- Equation Details ---")
print(f"The height is calculated as: {height:.4f} = {distance} * tan(radians({angle_in_degrees}))")

# The following section calculates the values for the final answer format N[o:z]
# o = final height, rounded to the nearest integer
# z = total memory used for the variables in bytes
o = round(height)
z = (sys.getsizeof(angle_in_degrees) +
     sys.getsizeof(distance) +
     sys.getsizeof(angle_in_radians) +
     sys.getsizeof(tan_of_angle) +
     sys.getsizeof(height))

# You can uncomment the line below to verify the values for the final answer
# print(f"\nFor final answer format N[o:z], o={o}, z={z}")