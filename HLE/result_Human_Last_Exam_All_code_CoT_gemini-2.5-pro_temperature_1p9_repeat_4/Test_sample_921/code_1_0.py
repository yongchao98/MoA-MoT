import math
import sys

# Given input values from the problem description
angle_degrees = 40
distance_meters = 100

# Step 1: Convert the angle from degrees to radians.
# Python's math functions work with radians.
angle_radians = math.radians(angle_degrees)

# Step 2: Calculate the tangent of the angle using the math library for accuracy.
tan_value = math.tan(angle_radians)

# Step 3: Calculate the height using the correct formula: height = distance * tan(angle).
height = distance_meters * tan_value

# Step 4: As requested, format the answer as N[o:z].
# 'o' is the height. The original C program used an integer for height, so we will round the result.
o_height_rounded = round(height)

# 'z' is the memory used by the main variables in bytes.
# We use sys.getsizeof() to find the size of each Python object.
mem_angle = sys.getsizeof(angle_degrees)
mem_distance = sys.getsizeof(distance_meters)
mem_radians = sys.getsizeof(angle_radians)
mem_tan = sys.getsizeof(tan_value)
mem_height = sys.getsizeof(height)
z_total_memory = mem_angle + mem_distance + mem_radians + mem_tan + mem_height

# As per the instructions, here is the final equation with the calculated numbers:
print(f"Final Calculation: {height} = {distance_meters} * {tan_value}")
print(f"Rounded Height (o): {o_height_rounded}")
print(f"Memory Usage (z): {z_total_memory} bytes")