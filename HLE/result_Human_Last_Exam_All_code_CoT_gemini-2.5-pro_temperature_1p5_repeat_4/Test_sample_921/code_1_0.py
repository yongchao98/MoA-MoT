import math
import sys

def calculate_palace_height():
    """
    Calculates the height of the palace accurately and prints the final equation.
    """
    # 1. Define input values from the problem description
    angle_deg = 40
    distance = 100

    # 2. Perform accurate calculation using the math library
    # Convert angle from degrees to radians for the tan function
    angle_rad = math.radians(angle_deg)

    # Calculate the tangent of the angle
    tan_val = math.tan(angle_rad)

    # Calculate the final height
    height = distance * tan_val

    # 3. As the original program output an integer, we round the result
    o = round(height)

    # 4. As requested, output each number in the final equation.
    # The equation is: height = distance * tan(angle)
    # We print the numbers involved in this final equation.
    print(f"The final equation with the calculated numbers is:")
    print(f"{o} = {distance} * tan({angle_deg})")

    # The variable 'z' for the final answer represents memory usage.
    # It is calculated below but not printed, as it is part of the final answer format.
    # The sizes are for a standard 64-bit Python environment.
    mem_angle_deg = sys.getsizeof(angle_deg)
    mem_distance = sys.getsizeof(distance)
    mem_angle_rad = sys.getsizeof(angle_rad)
    mem_tan_val = sys.getsizeof(tan_val)
    mem_height = sys.getsizeof(height)
    mem_o = sys.getsizeof(o)
    z = mem_angle_deg + mem_distance + mem_angle_rad + mem_tan_val + mem_height + mem_o

if __name__ == '__main__':
    calculate_palace_height()
