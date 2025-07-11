import math
import sys

def calculate_palace_height():
    """
    Calculates the height of a palace given the angle and distance,
    and also computes the memory usage of the variables.
    """
    # Given inputs
    angle_deg = 40
    distance = 100

    # To calculate the height accurately, we use floating-point arithmetic
    # and the standard math library, which is highly optimized.

    # 1. Convert the angle from degrees to radians
    angle_rad = math.radians(angle_deg)

    # 2. Calculate the tangent of the angle
    tangent_val = math.tan(angle_rad)

    # 3. Calculate the final height using the formula: height = distance * tan(angle)
    height = distance * tangent_val

    # 4. Calculate the memory used by all variables involved in the calculation
    # In Python, sys.getsizeof() returns the size of an object in bytes.
    mem_angle_deg = sys.getsizeof(angle_deg)
    mem_distance = sys.getsizeof(distance)
    mem_angle_rad = sys.getsizeof(angle_rad)
    mem_tangent_val = sys.getsizeof(tangent_val)
    mem_height = sys.getsizeof(height)
    total_mem = mem_angle_deg + mem_distance + mem_angle_rad + mem_tangent_val + mem_height

    # Print the explanation and the final equation with all numbers
    print("An accurate calculation using Python's math library:")
    print("-" * 50)
    print(f"The equation for the height is: height = distance * tan(angle_in_radians)")
    print("Plugging in the numbers:")
    # We show the numbers for each part of the final equation
    print(f"{height} = {distance} * {tangent_val}")
    print("-" * 50)

    # Output the final results for the answer format N[o:z]
    # 'o' is the calculated height (rounded to the nearest integer)
    # 'z' is the total memory used by the variables
    o = round(height)
    z = total_mem
    print(f"Calculated height (o) for the final answer: {o} meters")
    print(f"Total memory for variables (z) for the final answer: {z} bytes")

if __name__ == '__main__':
    calculate_palace_height()
