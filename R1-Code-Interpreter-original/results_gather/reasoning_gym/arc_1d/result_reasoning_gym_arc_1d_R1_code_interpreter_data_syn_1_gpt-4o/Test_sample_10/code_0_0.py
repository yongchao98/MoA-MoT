def transform_grid(input_grid):
    # Extract non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    # Calculate the number of zeros to append
    num_zeros = len(input_grid) - len(non_zero_elements)
    # Create the output grid
    output_grid = non_zero_elements + [0] * num_zeros
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 1, 5, 0, 0, 0]
# Transform the grid
output_grid = transform_grid(input_grid)
print(output_grid)