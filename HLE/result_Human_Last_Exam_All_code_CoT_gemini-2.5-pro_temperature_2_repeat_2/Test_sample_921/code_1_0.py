import math
import sys

def calculate_palace_height():
    """
    Calculates the height of the palace accurately and determines the memory
    usage of the variables involved.
    """
    # Input values
    angle_deg = 40
    distance_m = 100

    # 1. Convert angle from degrees to radians for trigonometric functions.
    angle_rad = math.radians(angle_deg)

    # 2. Calculate the tangent of the angle.
    tangent_val = math.tan(angle_rad)

    # 3. Calculate the height.
    # The final equation is: height = distance * tan(angle)
    height_m = distance_m * tangent_val

    # 4. Calculate the memory size for each variable used in the calculation.
    mem_angle_deg = sys.getsizeof(angle_deg)
    mem_distance_m = sys.getsizeof(distance_m)
    mem_angle_rad = sys.getsizeof(angle_rad)
    mem_tangent_val = sys.getsizeof(tangent_val)
    mem_height_m = sys.getsizeof(height_m)
    total_memory_bytes = mem_angle_deg + mem_distance_m + mem_angle_rad + mem_tangent_val + mem_height_m

    # 5. Format the output as N[o:z] and print it.
    # o = height, rounded to one decimal place.
    # z = total memory used in bytes.
    # The numbers in the final equation are: height_m, distance_m, and tangent_val.
    answer_string = f"N[{height_m:.1f}:{total_memory_bytes}]"
    print(answer_string)

calculate_palace_height()