import math
import sys

# Define the inputs given in the problem
angle_degrees = 40
distance_meters = 100

# Convert the angle from degrees to radians for the tan function
angle_radians = math.radians(angle_degrees)

# Calculate the tangent of the angle
tangent_value = math.tan(angle_radians)

# Calculate the height of the palace
height = distance_meters * tangent_value

# As requested, output each number in the final equation:
# height = distance * tangent_value
print(f"Calculated Height: {height}")
print(f"Initial Distance: {distance_meters}")
print(f"Calculated Tangent: {tangent_value}")

# Calculate the memory usage for all variables as per the request
# The variables are: angle_degrees, distance_meters, angle_radians, tangent_value, height
memory_angle = sys.getsizeof(angle_degrees)
memory_distance = sys.getsizeof(distance_meters)
memory_radians = sys.getsizeof(angle_radians)
memory_tangent = sys.getsizeof(tangent_value)
memory_height = sys.getsizeof(height)

total_memory = memory_angle + memory_distance + memory_radians + memory_tangent + memory_height

# The final answer needs to be constructed in the format N[o:z]
# o is the calculated height, rounded to the nearest integer
o = round(height)
# z is the total memory used for the variables
z = total_memory

# This final print generates the required answer string.
# Note: The problem asks me to return this in a special format,
# which I will do as the final line of my entire response.
# print(f"N[{o}:{z}]")