import math
import sys

# 1. Define input variables
angle_degrees = 40
distance = 100

# 2. Perform the calculation accurately
# Convert angle from degrees to radians
angle_radians = math.radians(angle_degrees)
# Calculate the tangent of the angle
tan_value = math.tan(angle_radians)
# Calculate the final height
height = distance * tan_value

# 3. Calculate memory usage
# We sum the memory (in bytes) of the primary variables used.
# The size can vary slightly based on the Python version and system architecture.
variables = [angle_degrees, distance, angle_radians, tan_value, height]
total_memory = sum(sys.getsizeof(v) for v in variables)

# 4. Output the results as requested
# The final equation requires printing each number
print(f"Corrected Calculation Steps:")
print(f"Height = distance * tan(angle_degrees)")
print(f"Height = {distance} * tan({angle_degrees})")
print(f"Height = {distance} * {tan_value}")
print(f"Final Height = {height} meters")

# Construct the final answer in the format N[o:z]
# o = calculated height
# z = total memory used for variables in bytes
final_answer_string = f"N[{height}:{total_memory}]"

print("\nFinal Answer String:")
print(final_answer_string)