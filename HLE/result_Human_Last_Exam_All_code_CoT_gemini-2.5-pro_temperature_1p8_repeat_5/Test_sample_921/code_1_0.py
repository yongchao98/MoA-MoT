import math
import sys

def calculate_palace_height():
    """
    Calculates the height of the palace accurately and determines the memory
    usage of the variables involved.
    """
    # Input values from the problem
    angle_deg = 40
    distance_m = 100

    # The accurate formula for the height is: height = distance * tan(angle)
    # 1. Convert angle from degrees to radians for the tan function
    angle_rad = math.radians(angle_deg)

    # 2. Calculate the tangent of the angle
    tangent_val = math.tan(angle_rad)

    # 3. Calculate the final height
    height_m = distance_m * tangent_val

    # As requested, output each number in the final equation:
    # height_m = distance_m * tangent_val
    print("Final Calculation:")
    print(f"{height_m} = {distance_m} * {tangent_val}")

    # Calculate the memory used for all variables in bytes
    mem_angle_deg = sys.getsizeof(angle_deg)
    mem_distance_m = sys.getsizeof(distance_m)
    mem_angle_rad = sys.getsizeof(angle_rad)
    mem_tangent_val = sys.getsizeof(tangent_val)
    mem_height_m = sys.getsizeof(height_m)
    total_memory = mem_angle_deg + mem_distance_m + mem_angle_rad + mem_tangent_val + mem_height_m

    # For the final answer format N[o:z]
    o = f"{height_m:.1f}"  # Calculated height, rounded for the answer string
    z = total_memory        # Total memory usage
    
    # This print statement is for generating the final answer content for the user.
    # It does not need to be part of the core logic, but helps present the solution.
    # print(f"Answer Format: N[{o}:{z}]")


if __name__ == '__main__':
    calculate_palace_height()
