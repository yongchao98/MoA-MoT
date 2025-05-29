def transform_grid(input_grid):
    # Find the start and end indices of the sequence of zeros
    zero_start = None
    zero_end = None
    for i in range(len(input_grid)):
        if input_grid[i] == 0:
            if zero_start is None:
                zero_start = i
            zero_end = i

    # Calculate the length of the zero sequence
    zero_length = zero_end - zero_start + 1

    # Shift the zero sequence to the left by 4 positions
    shift_amount = 4
    new_zero_start = max(0, zero_start - shift_amount)
    new_zero_end = new_zero_start + zero_length - 1

    # Create the output grid
    output_grid = input_grid[:]

    # Move the zeros
    for i in range(zero_length):
        output_grid[new_zero_start + i] = 0

    # Fill the remaining positions with the non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    non_zero_index = 0
    for i in range(len(output_grid)):
        if i < new_zero_start or i > new_zero_end:
            output_grid[i] = non_zero_elements[non_zero_index]
            non_zero_index += 1

    return output_grid

# Test the function with the given input
input_grid = [0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")