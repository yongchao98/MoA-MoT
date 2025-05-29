def transform_grid(input_grid):
    # Count the number of zeros
    zero_count = input_grid.count(0)
    # Count the number of non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    # Create the output grid
    output_grid = [0] * zero_count + non_zero_elements
    return output_grid

# Test input
input_grid = [6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)