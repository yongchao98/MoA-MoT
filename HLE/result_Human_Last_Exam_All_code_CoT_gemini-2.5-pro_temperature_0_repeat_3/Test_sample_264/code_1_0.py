from fractions import Fraction

def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two cubes shifted along the main diagonal.
    """
    # Define the properties of the cubes
    side_length = Fraction(1)
    shift = Fraction(1, 2)

    # Define the boundaries of the first cube (Cube A)
    cube_a_min = Fraction(0)
    cube_a_max = side_length

    # Define the boundaries of the second cube (Cube B)
    cube_b_min = shift
    cube_b_max = side_length + shift

    # Calculate the boundaries of the intersection for one dimension
    # The result is the same for all three dimensions (x, y, z)
    intersect_min = max(cube_a_min, cube_b_min)
    intersect_max = min(cube_a_max, cube_b_max)

    # Calculate the side length of the intersection cube
    intersect_side_length = intersect_max - intersect_min

    # Calculate the volume of the intersection
    volume = intersect_side_length ** 3

    # Print the final equation
    print(f"The intersection is a cube with side length {intersect_side_length}.")
    print("The volume is calculated as:")
    print(f"{intersect_side_length} * {intersect_side_length} * {intersect_side_length} = {volume}")
    # Also printing in decimal for clarity
    print(f"{float(intersect_side_length)} * {float(intersect_side_length)} * {float(intersect_side_length)} = {float(volume)}")


solve_cube_intersection()