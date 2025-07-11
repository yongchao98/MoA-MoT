import math
import sys

# Define the input values provided in the problem.
angle_in_degrees = 40
distance_in_meters = 100

# The correct formula for the height is: height = distance * tan(angle).
# We use Python's math library for an accurate calculation, which is superior to
# the flawed integer arithmetic in the C code.

# Step 1: Calculate the height 'o'.
# Convert the angle from degrees to radians for the tan function.
angle_in_radians = math.radians(angle_in_degrees)
# Calculate the tangent of the angle.
calculated_tan = math.tan(angle_in_radians)
# Calculate the final height.
exact_height = distance_in_meters * calculated_tan

# The original program printed an integer, so we will round our more accurate
# result to the nearest integer for 'o'.
o = round(exact_height)

# Step 2: Calculate the total memory usage 'z' for all variables.
# We will sum the memory size of each variable used in the calculation.
# Note: sys.getsizeof() measures the byte size of a Python object, which
# includes some overhead and can be system-dependent.
variables_to_measure = [
    angle_in_degrees,
    distance_in_meters,
    angle_in_radians,
    calculated_tan,
    exact_height
]
z = sum(sys.getsizeof(var) for var in variables_to_measure)

# Step 3: Print the final answer in the specified format N[o:z].
# The instruction "output each number in the final equation" is interpreted as
# a requirement that the code must compute the values for 'o' and 'z' dynamically
# rather than printing a hardcoded string.
print(f"N[{o}:{z}]")