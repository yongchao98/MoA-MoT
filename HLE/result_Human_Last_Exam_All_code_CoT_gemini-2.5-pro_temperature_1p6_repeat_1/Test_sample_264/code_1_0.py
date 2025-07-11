def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two unit cubes shifted
    by 1/2 along the main diagonal relative to each other.
    """
    # Side length of the two identical cubes
    side_length_L = 1
    
    # The shift is 1/2 along the main diagonal. We interpret this as a shift
    # vector of (1/2, 1/2, 1/2). So the shift along each axis is 1/2.
    shift_s = 0.5
    
    print(f"Let's assume the first cube occupies the region from 0 to {side_length_L} on each axis.")
    print(f"The second cube is shifted by {shift_s} along each axis.")
    print("-" * 30)

    # The intersection of the two cubes forms a smaller cube.
    # We calculate the side length of this intersection cube.
    # The interval for the first cube on an axis is [0, L].
    # The interval for the second cube on an axis is [s, L+s].
    # The intersection of [0, L] and [s, L+s] is [s, L].
    # The length of this intersection is L - s.
    
    intersect_side_length = side_length_L - shift_s
    
    print(f"The side length of the original cubes is L = {side_length_L}.")
    print(f"The shift distance along each axis is s = {shift_s}.")
    print(f"The side length of the intersection cube is L_intersect = L - s = {side_length_L} - {shift_s} = {intersect_side_length}.")
    print("-" * 30)

    # The volume of the intersection cube is its side length cubed.
    volume = intersect_side_length ** 3
    
    print("The volume of the intersection is the cube of its side length.")
    # The final equation is printed with each number, as requested.
    print(f"Volume = {intersect_side_length} * {intersect_side_length} * {intersect_side_length} = {volume}")

solve_cube_intersection()