import math
import sys

# The task is to determine the height given a fixed angle and distance.
# Angle = 40 degrees
# Distance = 100 meters

# Plan:
# 1. Use the correct trigonometric formula: height = distance * tan(angle).
# 2. Use the `math` library for accurate calculations of radians and tangent.
# 3. Calculate the total memory used by the variables involved in the calculation.
# 4. Print all the numbers involved as requested by the prompt.

# Given input values
angle_deg = 40
distance_m = 100

# 1. Calculate the accurate height
# First, convert the angle from degrees to radians for the tan function.
angle_rad = math.radians(angle_deg)
# Next, calculate the tangent of the angle.
tangent_val = math.tan(angle_rad)
# Finally, calculate the height.
height_float = distance_m * tangent_val
# The C program uses an integer for height, so we will round to the nearest whole number.
# This value corresponds to 'o' in the answer format N[o:z].
height_rounded = round(height_float)

# 2. Calculate the memory usage
# We use sys.getsizeof() to get the size of each Python object in bytes.
# We then sum the sizes for all variables used in our calculation.
mem_angle = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance_m)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tangent = sys.getsizeof(tangent_val)
mem_height = sys.getsizeof(height_float)

# This value corresponds to 'z' in the answer format N[o:z].
total_memory_bytes = mem_angle + mem_distance + mem_angle_rad + mem_tangent + mem_height

# 3. Output the required numbers from the final calculation.
# The "equation" is height = distance * tan(angle)
print(f"Distance (meters): {distance_m}")
print(f"Angle (degrees): {angle_deg}")
# Outputting 'o', the calculated height
print(f"Calculated Height 'o' (meters): {height_rounded}")
# Outputting 'z', the total memory used
print(f"Total Memory 'z' (bytes): {total_memory_bytes}")