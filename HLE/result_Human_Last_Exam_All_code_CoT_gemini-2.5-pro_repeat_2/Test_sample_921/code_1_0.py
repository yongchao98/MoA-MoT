import math
import sys

def calculate_palace_height():
    """
    Calculates the height of the palace accurately and determines the memory
    used for the variables.
    """
    # Inputs
    angle_deg = 40
    distance = 100

    # 1. Convert angle from degrees to radians for the math function
    angle_rad = math.radians(angle_deg)

    # 2. Calculate the tangent of the angle
    tangent_val = math.tan(angle_rad)

    # 3. Calculate the height
    height = distance * tangent_val
    
    # The accurate height is 'o' in the final answer, rounded to the nearest integer.
    o = round(height)

    # 4. Calculate memory usage for all variables
    # The variables are: angle_deg, distance, angle_rad, tangent_val, height
    mem_angle_deg = sys.getsizeof(angle_deg)
    mem_distance = sys.getsizeof(distance)
    mem_angle_rad = sys.getsizeof(angle_rad)
    mem_tangent_val = sys.getsizeof(tangent_val)
    mem_height = sys.getsizeof(height)
    
    # The total memory is 'z' in the final answer.
    z = mem_angle_deg + mem_distance + mem_angle_rad + mem_tangent_val + mem_height

    # As requested, output the numbers in the final equation
    print(f"Calculation: height({o}) = distance({distance}) * tan(angle_degrees({angle_deg}))")
    print(f"Total memory used for variables: {z} bytes")

    # 5. Print the final answer in the required format N[o:z]
    print(f"N[{o}:{z}]")

calculate_palace_height()