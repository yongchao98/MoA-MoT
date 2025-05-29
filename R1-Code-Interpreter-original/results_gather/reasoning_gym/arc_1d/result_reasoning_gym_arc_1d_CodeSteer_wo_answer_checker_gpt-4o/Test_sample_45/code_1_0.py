def transform_grid(input_grid):
    # Extract non-zero elements
    non_zero_elements = [num for num in input_grid if num != 0]
    # Calculate the number of zeros needed to fill the grid
    num_zeros = len(input_grid) - len(non_zero_elements)
    # Create the output grid
    output_grid = non_zero_elements + [0] * num_zeros
    return output_grid

# Test input
input_grid = [0, 8, 7, 0, 3, 4, 0, 3, 1, 0, 9, 0, 0, 0, 0, 9, 2, 7, 5, 4]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")