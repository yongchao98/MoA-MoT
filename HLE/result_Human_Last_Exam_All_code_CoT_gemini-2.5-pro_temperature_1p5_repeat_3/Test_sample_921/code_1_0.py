import math
import sys

def calculate_palace_height():
    """
    Calculates the height of the palace accurately and determines the memory
    used for the variables.
    """
    # Input values
    angle_degrees = 40
    distance = 100

    # For accurate calculations, we use floating-point numbers and the math library.
    # 1. Convert angle from degrees to radians.
    angle_radians = math.radians(angle_degrees)

    # 2. Calculate the height using the formula: height = distance * tan(angle).
    height = distance * math.tan(angle_radians)

    # 3. Round the height to the nearest integer for the output value 'o'.
    o = round(height)

    # 4. Calculate the memory used by all variables ('z').
    # We use sys.getsizeof() to find the size of each object in bytes.
    z = (sys.getsizeof(angle_degrees) +
         sys.getsizeof(distance) +
         sys.getsizeof(angle_radians) +
         sys.getsizeof(height))

    # Per the instructions, output each number in the final equation.
    # The equation is height = distance * tan(angle).
    print(f"{o} = {distance} * tan({angle_degrees})")
    
    # This print is for your information to see the final answer parts.
    # print(f"The answer format is N[o:z]")
    # print(f"Calculated o (rounded height): {o}")
    # print(f"Calculated z (memory in bytes): {z}")


calculate_palace_height()