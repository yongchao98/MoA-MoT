def find_square_combination_for_area():
    """
    Finds a combination of squares from the set S={2x2, 3x3, 5x5, 7x7}
    that can tile a target area. The specific target is the area of the
    smallest integer-length rectangle admitting a non-guillotine tiling,
    which is known to be a 12x12 square.
    """
    target_area = 144
    s7_area = 7 * 7
    s5_area = 5 * 5
    s3_area = 3 * 3
    s2_area = 2 * 2

    # Iterate through the number of largest squares to find a solution
    max_n7 = target_area // s7_area
    for n7 in range(max_n7, -1, -1):
        remaining_after_s7 = target_area - n7 * s7_area
        max_n5 = remaining_after_s7 // s5_area
        for n5 in range(max_n5, -1, -1):
            remaining_after_s5 = remaining_after_s7 - n5 * s5_area
            max_n3 = remaining_after_s5 // s3_area
            for n3 in range(max_n3, -1, -1):
                remaining_after_s3 = remaining_after_s5 - n3 * s3_area
                if remaining_after_s3 % s2_area == 0:
                    n2 = remaining_after_s3 // s2_area
                    
                    # Found a valid combination, print the equation
                    print(f"The area of the rectangle is {target_area}.")
                    print("One possible combination of squares that sums to this area is:")
                    print(f"({n7} * 7x7) + ({n5} * 5x5) + ({n3} * 3x3) + ({n2} * 2x2)")
                    
                    # Print the full equation with values as requested
                    print("\nThe area equation is:")
                    print(f"{target_area} = {n7} * {s7_area} + {n5} * {s5_area} + {n3} * {s3_area} + {n2} * {s2_area}")
                    return

find_square_combination_for_area()