import math
import sys

# Define the inputs provided in the problem
angle_deg = 40
distance = 100

# 1. Convert angle from degrees to radians for the math library
angle_rad = math.radians(angle_deg)

# 2. Calculate the tangent of the angle
tan_value = math.tan(angle_rad)

# 3. Calculate the height using the correct formula
height = distance * tan_value

# 4. Calculate the memory usage for the essential variables
mem_angle_deg = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tan_value = sys.getsizeof(tan_value)
mem_height = sys.getsizeof(height)
total_memory_bytes = mem_angle_deg + mem_distance + mem_angle_rad + mem_tan_value + mem_height

# Output each number in the final equation as requested
print("The accurate calculation is: height = distance * tan(angle)")
print(f"height = {distance} * tan({angle_deg})")
print(f"height = {distance} * {tan_value}")
print(f"Final height = {height} meters")

# Prepare the values for the final answer format N[o:z]
# 'o' is the height truncated to an integer
o = int(height)
# 'z' is the total memory used in bytes
z = total_memory_bytes

print("\n---")
print(f"For the final answer format N[o:z]:")
print(f"o (integer height) = {o}")
print(f"z (total memory for variables in bytes) = {z}")