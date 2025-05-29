def transform_grid(input_grid):
    # Step 1: Identify the non-zero sequence
    first_non_zero = None
    last_non_zero = None
    for i, value in enumerate(input_grid):
        if value != 0:
            if first_non_zero is None:
                first_non_zero = i
            last_non_zero = i

    # Step 2: Calculate the shift amount
    leading_zeros = first_non_zero
    non_zero_length = last_non_zero - first_non_zero + 1
    shift_amount = leading_zeros  # Based on the pattern observed

    # Step 3: Construct the output grid
    output_grid = [0] * len(input_grid)
    for i in range(non_zero_length):
        output_grid[shift_amount + i] = input_grid[first_non_zero + i]

    # Step 4: Output the result
    print("<<<" + " ".join(map(str, output_grid)) + ">>>")

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
transform_grid(input_grid)