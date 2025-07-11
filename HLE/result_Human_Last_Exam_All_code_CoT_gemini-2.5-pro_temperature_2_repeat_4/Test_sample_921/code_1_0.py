import math
import sys

# Plan:
# 1. Use the correct formula: height = distance * tan(angle).
# 2. Use Python's `math` library for accurate trigonometric functions.
# 3. Calculate the height for the given angle (40 degrees) and distance (100m).
# 4. Round the result to the nearest integer.
# 5. Calculate the memory used by the essential variables.
# 6. Print the numbers involved in the final equation as requested.

# Input values from the problem
angle_deg = 40
distance_m = 100

# Perform the accurate calculation
# The math.tan() function requires the angle to be in radians
height_float = distance_m * math.tan(math.radians(angle_deg))

# The original C code intended to produce an integer height.
# We round the accurate float result to the nearest whole number.
final_height_m = int(round(height_float))

# The problem requires printing each number in the final equation:
# final_height = distance * tan(angle)
print("Equation: height = distance * tan(angle)")
print(f"Calculated height: {final_height_m}")
print(f"Distance: {distance_m}")
print(f"Angle: {angle_deg}")

# Bonus: Calculate memory for the required format N[o:z]
# We calculate the size of all variables assigned in this script.
# On a typical 64-bit system, an int uses 28 bytes and a float uses 24 bytes.
mem_angle = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance_m)
mem_height_float = sys.getsizeof(height_float)
mem_final_height = sys.getsizeof(final_height_m)

total_memory_bytes = mem_angle + mem_distance + mem_height_float + mem_final_height

# This print is for context, the final answer will be provided below.
# print(f"Optimal Height (o): {final_height_m}")
# print(f"Memory Used (z): {total_memory_bytes} bytes")