import math
import sys

def calculate_palace_height():
    """
    Calculates the height of the palace accurately and determines the memory
    usage of the variables involved.

    The original C program is incorrect due to mathematical errors in the radian
    conversion, a flawed tangent approximation, and a critical integer overflow
    in the final calculation. This program provides the correct calculation.
    """
    # Input values from the problem description
    angle_degrees = 40
    distance = 100

    # Convert angle from degrees to radians for the tan function
    angle_radians = math.radians(angle_degrees)

    # Calculate the height using the correct trigonometric formula
    # height = distance * tan(angle)
    height_accurate = distance * math.tan(angle_radians)

    # 'o' is the calculated height, rounded to the nearest integer
    o = int(round(height_accurate))

    # 'z' is the total memory used for the primary variables in bytes.
    # sys.getsizeof() returns the base memory consumption of a Python object.
    mem_angle = sys.getsizeof(angle_degrees)
    mem_distance = sys.getsizeof(distance)
    mem_radians = sys.getsizeof(angle_radians)
    mem_height = sys.getsizeof(height_accurate)
    z = mem_angle + mem_distance + mem_radians + mem_height

    # The final equation is represented by the formatted string N[o:z].
    # The numbers in this equation are 'o' (the height) and 'z' (the memory).
    # The code prints the final answer in the required format.
    print(f"N[{o}:{z}]")

calculate_palace_height()