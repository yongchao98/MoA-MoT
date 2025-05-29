def transform_grid(input_grid):
    # Find the first and last non-zero indices
    first_non_zero = next(i for i, x in enumerate(input_grid) if x != 0)
    last_non_zero = len(input_grid) - next(i for i, x in enumerate(reversed(input_grid)) if x != 0) - 1
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero:last_non_zero + 1]
    
    # Calculate the new leading and trailing zeros
    total_length = len(input_grid)
    non_zero_length = len(non_zero_sequence)
    leading_zeros = (total_length - non_zero_length) // 2
    trailing_zeros = total_length - non_zero_length - leading_zeros
    
    # Construct the output grid
    output_grid = [0] * leading_zeros + non_zero_sequence + [0] * trailing_zeros
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 9, 9, 9, 9, 9, 9, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)