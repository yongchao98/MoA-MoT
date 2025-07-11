import math
import sys

# The original C program is incorrect due to mathematical errors and potential integer overflows.
# This Python script provides an accurate calculation.

# --- Inputs ---
# The given values from the problem description.
angle_deg = 40
distance = 100

# --- Calculation ---
# The correct formula is: height = distance * tan(angle)

# 1. Convert the angle from degrees to radians to use with the math library.
angle_rad = math.radians(angle_deg)

# 2. Calculate the tangent of the angle.
tan_value = math.tan(angle_rad)

# 3. Calculate the final height.
height = distance * tan_value

# --- Format the Answer ---
# The required format is N[o:z]
# 'o' is the calculated height, rounded to the nearest integer.
o = round(height)

# 'z' is the memory used for all the main variables, in bytes.
# Note: The result of sys.getsizeof() can vary slightly depending on the
# Python version and system architecture (e.g., 32-bit vs 64-bit).
# This is calculated on a standard 64-bit system.
variables = [angle_deg, distance, angle_rad, tan_value, height]
z = sum(sys.getsizeof(v) for v in variables)

# --- Output ---
# As requested, we output the numbers in the final equation.
print(f"The equation to solve is: height = distance * tan(angle_in_degrees)")
print(f"The values are: height = {distance} * tan({angle_deg})")
print(f"The result is: height = {height:.4f} meters")
print("-" * 20)
print(f"Calculated rounded height (o): {o}")
print(f"Memory for variables (z): {z} bytes")
print("-" * 20)
# The final answer is N[o:z] because the original C program is incorrect.
print(f"Final Answer String: N[{o}:{z}]")
