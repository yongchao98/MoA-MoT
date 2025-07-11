import math
import sys

# Define the input values
angle_deg = 40
distance = 100

# Convert the angle from degrees to radians for the calculation
angle_rad = math.radians(angle_deg)

# Calculate the tangent of the angle
tan_value = math.tan(angle_rad)

# Calculate the precise height of the palace
height = distance * tan_value

# Round the height to the nearest integer for the final result
final_height = int(round(height))

# Calculate the memory used for all the variables in bytes.
# These values are typical for a 64-bit Python implementation.
mem_angle = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tan = sys.getsizeof(tan_value)
mem_height = sys.getsizeof(height)

total_memory = mem_angle + mem_distance + mem_angle_rad + mem_tan + mem_height

# The final answer format is N[o:z].
# 'o' is the calculated height (final_height).
# 'z' is the total memory used for the variables (total_memory).
# We will present this in the required format at the very end.

# As requested, output each number in the final equation.
# The calculation is distance * tan(angle_degrees) = final_height
print(f"{distance} * tan({angle_deg}) = {final_height}")

# For clarity, here are the calculated values for the answer N[o:z]
# o = final_height
# z = total_memory
# print(f"Answer components: o = {final_height}, z = {total_memory}")