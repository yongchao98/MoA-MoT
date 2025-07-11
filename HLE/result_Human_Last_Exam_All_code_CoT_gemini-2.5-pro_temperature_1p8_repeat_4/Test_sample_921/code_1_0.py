import math
import sys

# Define the input values
angle_deg = 40
distance = 100

# 1. Convert angle from degrees to radians for the math functions
angle_rad = math.radians(angle_deg)

# 2. Calculate the tangent of the angle
tan_val = math.tan(angle_rad)

# 3. Calculate the height
height = distance * tan_val

# 4. As requested, output each number in the final equation
# The equation is: height = distance * tan(angle_in_radians)
print(f"Final Calculation: {height} = {distance} * {tan_val}")

# 5. Calculate the memory used for all variables
mem_angle = sys.getsizeof(angle_deg)
mem_dist = sys.getsizeof(distance)
mem_rad = sys.getsizeof(angle_rad)
mem_tan = sys.getsizeof(tan_val)
mem_height = sys.getsizeof(height)

total_memory = mem_angle + mem_dist + mem_rad + mem_tan + mem_height

# Prepare the components for the final answer N[o:z]
# o: the calculated height, rounded to the nearest integer
# z: the total memory used for all variables in bytes
optimal_height_int = round(height)
answer_string = f"N[{optimal_height_int}:{total_memory}]"

# This print statement is not part of the primary output but shows the answer
# format required by the user in the final response.
# print(f"Formatted Answer: {answer_string}")
