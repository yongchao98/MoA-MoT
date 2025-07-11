import math
import sys

# Define the input values from the problem
angle_deg = 40
distance = 100

# Convert the angle from degrees to radians for the tan function
angle_rad = math.radians(angle_deg)

# Calculate the tangent of the angle
tangent_val = math.tan(angle_rad)

# Calculate the final height
# height = distance * tan(angle)
height = distance * tangent_val

# As requested, output each number in the final equation
print(f"The correct calculation is:")
print(f"{distance} * tan({angle_deg}) = {height}")

# For the final answer string, 'o' is the height rounded to the nearest integer
o = round(height)

# 'z' is the total memory in bytes used for the main variables.
# Note: The size of objects in Python can vary slightly depending on the system.
# This calculation is based on a standard 64-bit system.
mem_angle_deg = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tangent_val = sys.getsizeof(tangent_val)
mem_height = sys.getsizeof(height)

z = mem_angle_deg + mem_distance + mem_angle_rad + mem_tangent_val + mem_height

# Print the final answer in the specified N[o:z] format
print("\nFinal Answer String:")
print(f"N[{o}:{z}]")