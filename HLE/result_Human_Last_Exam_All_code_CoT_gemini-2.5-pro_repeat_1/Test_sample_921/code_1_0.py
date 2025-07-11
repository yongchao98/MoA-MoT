import math
import sys

def calculate_palace_height():
    """
    Calculates the height of the palace accurately and determines the memory
    usage of the variables involved.
    """
    # Input values from the problem
    angle_deg = 40
    distance = 100

    # For the most accurate calculation, we use floating-point arithmetic.
    # 1. Convert angle from degrees to radians
    angle_rad = math.radians(angle_deg)

    # 2. Calculate the tangent of the angle
    tan_value = math.tan(angle_rad)

    # 3. Calculate the height
    height = distance * tan_value

    # As requested, output each number in the final equation.
    # The equation is: height = distance * tan(angle)
    # The numbers are the distance, the angle, and the final calculated height.
    print(f"{distance} * tan({angle_deg}) = {round(height)}")

    # The final answer N[o:z] requires the memory usage 'z'.
    # We calculate it here for the final response, but don't print it.
    # z = sum of memory for all variables used in the calculation.
    mem_angle_deg = sys.getsizeof(angle_deg)
    mem_distance = sys.getsizeof(distance)
    mem_angle_rad = sys.getsizeof(angle_rad)
    mem_tan_value = sys.getsizeof(tan_value)
    mem_height = sys.getsizeof(height)
    
    total_memory = mem_angle_deg + mem_distance + mem_angle_rad + mem_tan_value + mem_height
    # o = round(height) = 84
    # z = total_memory = 28 + 28 + 24 + 24 + 24 = 128 (on a typical 64-bit system)

if __name__ == '__main__':
    calculate_palace_height()
