import math
import sys

# Given values from the problem
angle_deg = 40
distance_m = 100

# Perform the calculation accurately using floating-point math
# 1. Convert angle from degrees to radians
angle_rad = math.radians(angle_deg)

# 2. Calculate the tangent of the angle
tan_of_angle = math.tan(angle_rad)

# 3. Calculate the height
height = distance_m * tan_of_angle

# The required answer format is N[o:z]
# o: The calculated height, rounded to the nearest integer.
o = round(height)

# z: The memory used for all the variables involved in the calculation.
# We sum the memory size of each variable using sys.getsizeof().
# This assumes a standard 64-bit Python interpreter.
mem_angle = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance_m)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tan = sys.getsizeof(tan_of_angle)
mem_height = sys.getsizeof(height)

z = mem_angle + mem_distance + mem_angle_rad + mem_tan + mem_height

# As requested, output each number in the final equation.
# This symbolically represents the calculation performed.
print(f"{o} = {distance_m} * tan({angle_deg})")

# Print the final result in the specified format 'N[o:z]'
print(f"N[{o}:{z}]")
