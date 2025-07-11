import math
import sys

# Define the input values
angle_deg = 40
distance = 100

# Convert angle from degrees to radians
angle_rad = math.radians(angle_deg)

# Calculate the tangent of the angle
tangent_val = math.tan(angle_rad)

# Calculate the height
height = distance * tangent_val

# Calculate the memory used by the primary variables in bytes
mem_angle_deg = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tangent_val = sys.getsizeof(tangent_val)
mem_height = sys.getsizeof(height)
total_memory = mem_angle_deg + mem_distance + mem_angle_rad + mem_tangent_val + mem_height

# Round the final height to the nearest integer
rounded_height = round(height)

# Output each number in the final equation as requested
print(f"Calculation: height = distance * tan(angle_deg)")
print(f"Values: height = {distance} * {tangent_val}")
print(f"Result: height = {height}")

# The final answer format is N[o:z] where o is the optimal height
# and z is the memory used for all variables.
# We print this result so you don't have to copy it.
print("\n---")
print("Final Answer:")
print(f"N[{rounded_height}:{total_memory}]")