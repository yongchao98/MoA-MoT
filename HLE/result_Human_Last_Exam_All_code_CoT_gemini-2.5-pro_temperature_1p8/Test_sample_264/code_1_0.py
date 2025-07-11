def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two cubes shifted along the main diagonal.
    """
    # Side length of the original cubes
    original_side = 1.0

    # The shift along each axis (x, y, and z)
    shift = 0.5

    # The intersection of the two cubes forms a new, smaller cube.
    # We calculate the side length of this intersection cube.
    # For any axis, the first cube is in the range [0, 1].
    # The second cube is in the range [0.5, 1.5].
    # The intersection is the range [0.5, 1], which has a length of 1 - 0.5.
    intersection_side_length = original_side - shift

    # The volume is the side length of the intersection cubed.
    volume = intersection_side_length ** 3

    print(f"The original side length of the cubes is {original_side}.")
    print(f"The shift along the main diagonal is equivalent to a shift of {shift} along each axis.")
    print(f"The side length of the resulting intersection cube is the original side length minus the shift.")
    print(f"Side length of intersection = {original_side} - {shift} = {intersection_side_length}")
    print("\nThe volume of the intersection is calculated as (side length)^3.")
    print(f"Volume = {intersection_side_length} * {intersection_side_length} * {intersection_side_length} = {volume}")

solve_cube_intersection()