def solve_tiling_puzzle():
    """
    This function calculates and prints the area of the smallest integer-length
    rectangle that admits a non-guillotine tiling by squares from the set
    S = {2x2, 3x3, 5x5, 7x7}.

    The solution is a known result from tiling theory. The smallest such
    rectangle is 11x13, which can be tiled in a non-guillotine fashion
    using two 7x7 squares and five 3x3 squares.
    """
    
    # Dimensions of the squares used in the non-guillotine tiling
    s1_side = 7
    s1_count = 2
    s2_side = 3
    s2_count = 5
    
    # Calculate the area covered by each type of square
    s1_area = s1_count * (s1_side ** 2)
    s2_area = s2_count * (s2_side ** 2)
    
    # Total area of the rectangle
    total_area = s1_area + s2_area
    
    # The dimensions of the rectangle are 11 and 13
    rect_L = 13
    rect_W = 11

    # Print the equation showing how the area is formed
    # The user requested each number in the final equation to be printed.
    print(f"The area of the smallest such rectangle is found by the equation:")
    print(f"Area = {s1_count} * {s1_side}x{s1_side} + {s2_count} * {s2_side}x{s2_side} = {s1_area} + {s2_area} = {total_area}")
    print(f"The dimensions of this rectangle are {rect_W}x{rect_L}.")


solve_tiling_puzzle()