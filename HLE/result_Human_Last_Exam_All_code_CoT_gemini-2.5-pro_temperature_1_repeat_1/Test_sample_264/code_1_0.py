import math

def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two unit cubes
    shifted by 1/2 along the main diagonal.
    """
    
    # Side length of the original cubes
    side_length = 1.0
    
    # The shift value along each axis, based on the problem description.
    # A shift of 1/2 along the main diagonal (1,1,1) is interpreted
    # as a shift vector of (1/2, 1/2, 1/2).
    shift = 0.5
    
    # Define the bounds of the first cube (C1)
    c1_min = 0.0
    c1_max = side_length
    
    # Define the bounds of the second, shifted cube (C2)
    c2_min = shift
    c2_max = side_length + shift
    
    # The intersection of two intervals [a, b] and [c, d] is [max(a, c), min(b, d)].
    # The length of the intersection is min(b, d) - max(a, c).
    # Since the setup is symmetrical, the side length of the intersection
    # is the same for x, y, and z axes.
    
    intersection_side_length = min(c1_max, c2_max) - max(c1_min, c2_min)
    
    # The volume of the intersection is the cube of its side length.
    volume = intersection_side_length ** 3
    
    print("This script calculates the volume of the intersection of two shifted unit cubes.")
    print(f"The side length of each original cube is {side_length}.")
    print(f"The cubes are shifted relative to each other by {shift} along each axis.")
    print("-" * 30)
    print(f"The side length of the resulting intersection cube is {intersection_side_length}.")
    print("The volume is calculated by cubing this side length.")
    print("\nFinal Equation:")
    # The final equation shows the numbers used for the calculation.
    print(f"({intersection_side_length}) ** 3 = {volume}")

solve_cube_intersection()
<<<0.125>>>