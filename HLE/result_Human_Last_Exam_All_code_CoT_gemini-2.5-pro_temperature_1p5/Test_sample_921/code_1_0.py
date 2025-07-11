import math
import sys

# The final equation is height = distance * tan(angle)
angle_deg = 40
distance_m = 100

# 1. Calculate the values for the equation
# Python's math library provides accurate trigonometric functions.
angle_rad = math.radians(angle_deg)
tan_of_angle = math.tan(angle_rad)
height = distance_m * tan_of_angle

# As requested, here are the numbers in the final equation:
# The calculation is height = distance * tan(angle_in_radians)
print(f"height ({height:.4f}) = distance ({distance_m}) * tan_of_angle ({tan_of_angle:.4f})")


# 2. Prepare the final answer string N[o:z]
# o is the height rounded to the nearest integer
o = round(height)

# z is the memory used for all variables involved in the calculation
# sys.getsizeof() returns the memory footprint of an object in bytes.
# Variables: angle_deg, distance_m, angle_rad, tan_of_angle, height
mem_angle_deg = sys.getsizeof(angle_deg)
mem_distance_m = sys.getsizeof(distance_m)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tan_of_angle = sys.getsizeof(tan_of_angle)
mem_height = sys.getsizeof(height)

z = mem_angle_deg + mem_distance_m + mem_angle_rad + mem_tan_of_angle + mem_height

# The original C program is incorrect. The final answer format should be N[o:z].
print(f"\nN[{o}:{z}]")