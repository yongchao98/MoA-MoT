def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two unit cubes
    shifted by 1/2 along the main diagonal.
    """
    # The side length of each original cube.
    side_length = 1.0

    # We interpret the shift "by 1/2 along the main diagonal" as a shift
    # of 1/2 along each of the x, y, and z axes.
    shift_along_axis = 0.5

    # The first cube (C1) can be defined in the region:
    # 0 <= x <= 1, 0 <= y <= 1, 0 <= z <= 1.
    
    # The second cube (C2) is shifted from C1 and occupies the region:
    # 0.5 <= x <= 1.5, 0.5 <= y <= 1.5, 0.5 <= z <= 1.5.

    # The intersection is a smaller cube. Its side length is determined
    # by the overlap of the cubes' ranges on any axis.
    # The overlap of [0, 1] and [0.5, 1.5] is [0.5, 1].
    # The length of this overlapping interval is 1 - 0.5 = 0.5.
    intersection_side_length = side_length - shift_along_axis

    # The volume of the intersection is the cube of its side length.
    volume = intersection_side_length ** 3

    print("Step 1: Define the problem.")
    print(f"Two cubes of side length {side_length} are shifted relative to each other.")
    print(f"The shift along each axis is {shift_along_axis}.\n")

    print("Step 2: Find the dimensions of the intersection.")
    print(f"The intersection forms a smaller cube. Its side length is the overlap along any axis.")
    print(f"Side length of intersection = {side_length} - {shift_along_axis} = {intersection_side_length}\n")

    print("Step 3: Calculate the volume of the intersection.")
    print("The volume is the side length of the intersection cubed.")
    # The final equation with each number explicitly shown
    print(f"Volume = {intersection_side_length} * {intersection_side_length} * {intersection_side_length} = {volume}")

solve_cube_intersection()