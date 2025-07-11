def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two cubes of side length 1,
    where one is shifted by 1/2 along the main diagonal relative to the other.
    """
    # Side length of the original cubes
    side_length = 1.0

    # We define the first cube (C1) from the origin.
    # Its coordinates range from 0 to 1 on each axis.
    c1_min = 0.0
    c1_max = side_length

    # The shift is 1/2 along the main diagonal. We interpret this as a shift
    # of 1/2 along each axis.
    shift = 0.5

    # The second cube (C2) is shifted from the first one.
    # Its coordinates range from 0.5 to 1.5 on each axis.
    c2_min = shift
    c2_max = side_length + shift

    # The intersection is a smaller cube. Its range on each axis is found by
    # taking the maximum of the minimums and the minimum of the maximums.
    intersect_min = max(c1_min, c2_min)
    intersect_max = min(c1_max, c2_max)

    # The side length of the intersection cube.
    intersect_side_length = intersect_max - intersect_min

    # The volume of the intersection cube.
    volume = intersect_side_length ** 3
    
    # Print the final equation with each number.
    print(f"The volume of the intersection is the result of the following calculation:")
    print(f"{intersect_side_length} * {intersect_side_length} * {intersect_side_length} = {volume}")

solve_cube_intersection()