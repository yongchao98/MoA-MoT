import math
import sys

def calculate_palace_height():
    """
    Calculates the height of a palace accurately and determines the memory
    usage of the variables involved.
    """
    # Input values from the problem
    angle_deg = 40
    distance = 100

    # 1. Convert angle from degrees to radians for the tan function
    angle_rad = math.radians(angle_deg)

    # 2. Calculate the tangent of the angle
    tan_val = math.tan(angle_rad)

    # 3. Calculate the final height
    # The equation is: height = distance * tan(angle_in_radians)
    height = distance * tan_val
    
    # Per the output format requirements, 'o' is the integer part of the height.
    o = int(height)

    # Calculate 'z', the total memory used by the primary variables.
    # Note: sys.getsizeof() reports memory usage in bytes, which can vary
    # slightly depending on the Python version and system architecture (e.g., 32-bit vs 64-bit).
    # This calculation is based on a typical 64-bit system.
    mem_angle_deg = sys.getsizeof(angle_deg)
    mem_distance = sys.getsizeof(distance)
    mem_angle_rad = sys.getsizeof(angle_rad)
    mem_tan_val = sys.getsizeof(tan_val)
    mem_height = sys.getsizeof(height)
    z = mem_angle_deg + mem_distance + mem_angle_rad + mem_tan_val + mem_height

    # Print the equation with the numbers used, as requested.
    print(f"The calculation is: height = distance * tan(angle_in_degrees)")
    print(f"Substituting the values: {o} = {distance} * tan({angle_deg})")
    
    # Final values for the answer format N[o:z]
    print(f"\nCalculated height 'o': {o}")
    print(f"Memory used for variables 'z': {z} bytes")

calculate_palace_height()