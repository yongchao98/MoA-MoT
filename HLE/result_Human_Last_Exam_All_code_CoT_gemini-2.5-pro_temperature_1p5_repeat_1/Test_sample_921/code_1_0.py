import math
import sys

def calculate_palace_height():
    """
    Calculates the height of the palace with high accuracy and reports memory usage.
    """
    # Input values
    angle_deg = 40
    distance_m = 100

    # Convert angle from degrees to radians
    angle_rad = math.radians(angle_deg)

    # Calculate the tangent of the angle
    tan_val = math.tan(angle_rad)

    # Calculate the height
    height_m = distance_m * tan_val

    # Calculate the memory used by all relevant variables in bytes
    # Note: Memory size can vary slightly depending on the Python version and system architecture (32-bit vs 64-bit).
    mem_angle = sys.getsizeof(angle_deg)
    mem_distance = sys.getsizeof(distance_m)
    mem_angle_rad = sys.getsizeof(angle_rad)
    mem_tan = sys.getsizeof(tan_val)
    mem_height = sys.getsizeof(height_m)
    total_memory = mem_angle + mem_distance + mem_angle_rad + mem_tan + mem_height

    # Per the instructions, output each number in the final equation.
    # The final equation is: height = distance * tan(angle_in_radians)
    print(f"{height_m} = {distance_m} * {tan_val}")

    # Output the final answer in the specified N[o:z] format
    print(f"N[{height_m}:{total_memory}]")

calculate_palace_height()