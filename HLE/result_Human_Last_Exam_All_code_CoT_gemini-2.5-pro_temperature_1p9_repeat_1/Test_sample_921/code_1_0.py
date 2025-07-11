import math
import sys

# This program provides an accurate calculation for the height of the palace.
# It also calculates the memory footprint of its variables as requested.

# --- Inputs ---
angle_degrees = 40
distance = 100

# --- Calculation ---

# Step 1: Convert the angle from degrees to radians.
angle_radians = math.radians(angle_degrees)

# Step 2: Calculate the tangent of the angle using the accurate math library function.
tangent_value = math.tan(angle_radians)

# Step 3: Calculate the height.
# The formula is: height = distance * tan(angle).
height = distance * tangent_value

# As requested, here are the numbers in the final equation:
print(f"height = distance * tan(angle_in_radians)")
print(f"{height} = {distance} * {tangent_value}")


# --- Final Answer Generation (for N[o:z] format) ---

# o: The height calculated by the program. We truncate it to an integer
#    to match the integer output format of the original C program.
o = int(height)

# z: The total memory used for all the main variables, in bytes.
mem_angle_degrees = sys.getsizeof(angle_degrees)
mem_distance = sys.getsizeof(distance)
mem_angle_radians = sys.getsizeof(angle_radians)
mem_tangent_value = sys.getsizeof(tangent_value)
mem_height = sys.getsizeof(height)
z = mem_angle_degrees + mem_distance + mem_angle_radians + mem_tangent_value + mem_height

# The constructed answer string is N[83:128] on a typical 64-bit system.
# This part is for analysis; the final answer is provided separately below.
# print(f"Answer Format: N[{o}:{z}]")
