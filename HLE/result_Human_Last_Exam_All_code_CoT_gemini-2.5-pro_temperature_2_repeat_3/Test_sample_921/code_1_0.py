import math
import sys

# The king wants to know the height of his highest palace.

# --- Inputs ---
# The angle and distance are given.
angle_deg = 40
distance = 100

# --- Accurate Calculation ---
# 1. Convert the angle from degrees to radians for the trigonometric function.
angle_rad = math.radians(angle_deg)

# 2. Calculate the tangent of the angle.
#    This uses a highly accurate floating-point implementation.
tan_of_angle = math.tan(angle_rad)

# 3. Calculate the final height.
#    height = distance * tan(angle)
height = distance * tan_of_angle

# --- Memory Usage Calculation ---
# Calculate the memory used by all variables in the calculation process.
# In Python, the size of an integer is typically 28 bytes and a float is 24 bytes.
size_angle_deg = sys.getsizeof(angle_deg)
size_distance = sys.getsizeof(distance)
size_angle_rad = sys.getsizeof(angle_rad)
size_tan_of_angle = sys.getsizeof(tan_of_angle)
size_height = sys.getsizeof(height)

total_memory = size_angle_deg + size_distance + size_angle_rad + size_tan_of_angle + size_height

# --- Output ---
# As requested, output each number in the final equation.
# The final equation is: height = distance * tan_of_angle
print(f"The final calculation is: {height} = {distance} * {tan_of_angle}")

# The final answer in the specified format N[o:z]
# 'o' is the height my program calculates (integer part).
# 'z' is the memory used for all my variables.
output_o = int(height)
output_z = total_memory

print(f"<<<N[{output_o}:{output_z}]>>>")