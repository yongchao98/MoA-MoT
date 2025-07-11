import math
import sys

# The user's C program is incorrect. We provide a correct solution here.
# This program calculates the height accurately and determines the memory used by its variables.

# Step 1: Define the given input values
angle_deg = 40
distance = 100

# Step 2: Perform the accurate calculation using the math library
# Convert angle from degrees to radians for the trigonometric function
angle_rad = math.radians(angle_deg)

# Calculate the tangent of the angle
tan_val = math.tan(angle_rad)

# The formula for the height is distance * tan(angle)
height = distance * tan_val

# Step 3: Calculate the total memory used for the variables
# Note: Memory size can vary depending on the system (e.g., 32-bit vs 64-bit).
# This is calculated for a typical 64-bit system.
mem_angle_deg = sys.getsizeof(angle_deg)
mem_distance = sys.getsizeof(distance)
mem_angle_rad = sys.getsizeof(angle_rad)
mem_tan_val = sys.getsizeof(tan_val)
mem_height = sys.getsizeof(height)
total_memory = mem_angle_deg + mem_distance + mem_angle_rad + mem_tan_val + mem_height

# Step 4: Output the results as requested by the prompt

# Output each number in the final equation. We show the result, and the two numbers multiplied.
print("Final Equation: height = distance * tan(radians(angle))")
print(f"{height:.2f} = {distance} * {tan_val:.2f}")

# Output the final answer string 'N[o:z]', where 'o' is the calculated
# height and 'z' is the total memory used by the variables in bytes.
print(f"N[{height:.2f}:{total_memory}]")