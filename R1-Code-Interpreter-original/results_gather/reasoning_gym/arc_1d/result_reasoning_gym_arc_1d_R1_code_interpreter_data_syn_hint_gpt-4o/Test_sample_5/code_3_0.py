def transform_grid(input_grid):
    # Find the first and last non-zero indices
    first_non_zero_index = next(i for i, x in enumerate(input_grid) if x != 0)
    last_non_zero_index = len(input_grid) - next(i for i, x in enumerate(reversed(input_grid)) if x != 0) - 1
    
    # Calculate the number of leading zeros in the output
    leading_zeros = first_non_zero_index - 3
    
    # Extract the non-zero segment
    non_zero_segment = input_grid[first_non_zero_index:last_non_zero_index + 1]
    
    # Calculate the number of trailing zeros in the output
    trailing_zeros = len(input_grid) - leading_zeros - len(non_zero_segment)
    
    # Construct the output grid
    output_grid = [0] * leading_zeros + non_zero_segment + [0] * trailing_zeros
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Transform the grid
output_grid = transform_grid(input_grid)
print(output_grid)