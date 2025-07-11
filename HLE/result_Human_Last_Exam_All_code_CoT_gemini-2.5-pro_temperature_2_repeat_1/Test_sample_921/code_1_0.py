import math
import sys

# Given input values
angle_deg = 40
distance = 100

# The formula to find the height is: height = distance * tan(angle)
# We need the angle in radians for the math.tan() function.

# Perform an accurate calculation using the math library
angle_rad = math.radians(angle_deg)
tan_val = math.tan(angle_rad)
height = distance * tan_val

# The problem requires outputting each number in the final equation.
# The rounded height 'o' will be used in the final equation display.
o = round(height)
print(f"Based on the accurate calculation, the final equation with the numbers is:")
print(f"{o} = {distance} * tan({angle_deg})")

# To determine 'z' for the final answer N[o:z], we sum the memory
# size of all the variables involved in the calculation.
z = (sys.getsizeof(angle_deg) +
     sys.getsizeof(distance) +
     sys.getsizeof(angle_rad) +
     sys.getsizeof(tan_val) +
     sys.getsizeof(height))

# The final answer needs to be printed in the N[o:z] format.
# Note: This is part of the program's output, representing the solution.
print("\nThe final answer string is:")
print(f"N[{o}:{z}]")