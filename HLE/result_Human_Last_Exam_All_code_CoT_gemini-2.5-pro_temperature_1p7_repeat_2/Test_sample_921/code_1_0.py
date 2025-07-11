import math
import sys

# Given values
angle_deg = 40
distance_m = 100

# 1. Convert angle from degrees to radians for the trigonometric function
angle_rad = math.radians(angle_deg)

# 2. Calculate the tangent of the angle
# The core formula for height is: height = distance * tan(angle)
tangent_val = math.tan(angle_rad)

# 3. Calculate the height
height = distance_m * tangent_val

# 4. Calculate the 'o' and 'z' values for the final answer format N[o:z]
# 'o' is the calculated height, rounded to the nearest integer
o = int(round(height))

# 'z' is the memory used for all primary variables in bytes
# Note: Memory size can vary slightly based on the Python version and OS.
z = sys.getsizeof(angle_deg) + sys.getsizeof(distance_m) + sys.getsizeof(angle_rad) + sys.getsizeof(tangent_val) + sys.getsizeof(height)

# Print the final equation with all its numbers as requested
print(f"The equation to find the height is: height = distance * tan(angle)")
print(f"Plugging in the numbers: height = {distance_m} * tan({angle_deg})")
print(f"Calculated height: {height:.2f} meters")
print(f"Rounded height (o): {o} meters")
print(f"Memory for variables (z): {z} bytes")
print("")
# Print the final answer in the required format
print(f"N[{o}:{z}]")
