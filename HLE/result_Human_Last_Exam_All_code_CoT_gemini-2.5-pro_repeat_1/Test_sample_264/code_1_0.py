def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two cubes of side length 1,
    shifted by 1/2 along the main diagonal relative to each other.
    """
    # Define the parameters from the problem description.
    side_length = 1.0
    shift = 0.5

    # The intersection of the two cubes forms a smaller cube.
    # Its side length is the original side length minus the shift along any one axis.
    intersection_side_length = side_length - shift

    # The volume of this intersection cube is its side length cubed.
    volume = intersection_side_length ** 3

    # Print the explanation and the final equation.
    print("The problem describes two cubes of side length 1, shifted by 1/2 along the main diagonal.")
    print("This creates an overlapping region which is also a cube.")
    print(f"The side length of this intersection cube is the original side length ({side_length}) minus the shift ({shift}).")
    print(f"Side length of intersection = {side_length} - {shift} = {intersection_side_length}")
    print("\nThe volume of the intersection is calculated by cubing its side length.")
    print(f"Volume = {intersection_side_length} * {intersection_side_length} * {intersection_side_length} = {volume}")

solve_cube_intersection()