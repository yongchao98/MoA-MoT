import math
import sys

# The original C program is incorrect. It contains errors in the radian conversion formula,
# the tangent approximation logic, and will suffer from a massive integer overflow during
# calculation. This Python program correctly calculates the palace height.

# --- Inputs ---
angle_deg = 40
distance = 100

print("Calculating height for a palace with a distance of 100m and an elevation angle of 40 degrees.")
print(f"The equation to solve is: height = distance * tan(angle)")
print(f"Substituted values: height = {distance} * tan({angle_deg} degrees)\n")

# --- Step-by-step Calculation ---

# Step 1: Convert the angle from degrees to radians for the tan function.
angle_rad = math.radians(angle_deg)
print(f"Step 1: Convert angle to radians -> {angle_rad}")

# Step 2: Calculate the tangent of the angle in radians.
tan_val = math.tan(angle_rad)
print(f"Step 2: Calculate tangent -> {tan_val}")

# Step 3: Calculate the final height.
# The final equation is height = 100 * tan(0.6981...)
# height = 100 * 0.8390...
height = distance * tan_val
print(f"Step 3: Calculate height -> {int(height)} = {distance} * {tan_val}")


# --- Final Answer Formatting ---

# Get the integer part of the height for the output 'o'.
output_height = int(height)

# Calculate the total memory in bytes used by the main variables for the output 'z'.
mem_angle_deg = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tan_val = sys.getsizeof(tan_val)
mem_height = sys.getsizeof(height)
total_memory = mem_angle_deg + mem_distance + mem_angle_rad + mem_tan_val + mem_height

# The C program is incorrect ('N'). 'o' is the height and 'z' is the memory used.
print(f"\nFinal Answer: N[{output_height}:{total_memory}]")