def transform_grid(input_grid):
    # Extract non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    # Calculate the number of zeros
    num_zeros = len(input_grid) - len(non_zero_elements)
    # Create the output grid by shifting non-zero elements to the right
    output_grid = [0] * num_zeros + non_zero_elements
    return output_grid

# Test input
input_grid = [6, 3, 4, 2, 6, 9, 8, 1, 7, 0, 0, 0]
# Transform the grid
output_grid = transform_grid(input_grid)
print(output_grid)