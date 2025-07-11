import math
import sys

def calculate_palace_height():
    """
    Calculates the palace height accurately and determines the memory usage
    of the variables involved.
    """
    # Input values from the problem
    angle_deg = 40
    distance = 100

    # 1. Convert angle from degrees to radians for the tan function
    angle_rad = math.radians(angle_deg)

    # 2. Calculate the height using the accurate math.tan function
    # The formula is: height = distance * tan(angle)
    accurate_height = distance * math.tan(angle_rad)
    
    # 3. Calculate the memory used by the key variables in bytes.
    # We use sys.getsizeof() to measure the size of the Python objects for
    # our inputs, the intermediate radian value, and the final result.
    mem_angle = sys.getsizeof(angle_deg)
    mem_distance = sys.getsizeof(distance)
    mem_rad = sys.getsizeof(angle_rad)
    mem_height = sys.getsizeof(accurate_height)

    total_memory = mem_angle + mem_distance + mem_rad + mem_height

    # 4. Format the final output string as N[o:z]
    # o = accurate height, rounded to two decimal places
    # z = total memory used for the variables
    print(f"N[{accurate_height:.2f}:{total_memory}]")

calculate_palace_height()