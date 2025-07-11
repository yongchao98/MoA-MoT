import math
import sys

# The final equation to find the height is: height = distance * tan(angle)
# The numbers for this equation are provided here.
angle_deg = 40
distance = 100

# To calculate the height accurately, we use Python's math library.
# First, convert the angle from degrees to radians.
angle_rad = math.radians(angle_deg)

# Second, calculate the tangent of the angle in radians.
tangent_val = math.tan(angle_rad)

# Finally, calculate the height.
height = distance * tangent_val

# The original C program is incorrect, so the answer starts with 'N'.
# 'o' is the accurate height, rounded to the nearest integer.
o = round(height)

# 'z' is the total memory in bytes used for the primary variables in this calculation.
# We use sys.getsizeof() to measure the memory footprint of each variable.
z = (sys.getsizeof(angle_deg) +
     sys.getsizeof(distance) +
     sys.getsizeof(angle_rad) +
     sys.getsizeof(tangent_val) +
     sys.getsizeof(height))

# The final output is printed in the specified format N[o:z].
print(f"N[{o}:{z}]")