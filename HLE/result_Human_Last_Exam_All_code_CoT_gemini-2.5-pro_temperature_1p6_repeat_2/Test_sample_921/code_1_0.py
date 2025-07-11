import math
import sys

# Define the inputs given in the problem
angle_deg = 40
distance = 100

# Convert the angle from degrees to radians for the trigonometric function
angle_rad = math.radians(angle_deg)

# Calculate the tangent of the angle
tangent_val = math.tan(angle_rad)

# Calculate the height of the palace
# The result is a float, so we round it to the nearest integer for the final height
height_float = distance * tangent_val
height_int = int(round(height_float))

# The problem asks to output each number in the final equation.
# The equation is: height = distance * tan(angle)
print("Calculation Breakdown:")
print(f"height = {distance} * tan({angle_deg})")
print(f"height = {height_int} meters (rounded)")

# Calculate the total memory used for all variables in the calculation
# The variables are: angle_deg, distance, angle_rad, tangent_val, height_float, height_int
variables = [angle_deg, distance, angle_rad, tangent_val, height_float, height_int]
total_memory_used = sum(sys.getsizeof(v) for v in variables)

# The requested format is N[o:z]
# o = calculated height
# z = memory used
print("\nFinal Answer Components:")
print(f"o (height): {height_int}")
print(f"z (memory): {total_memory_used} bytes")
