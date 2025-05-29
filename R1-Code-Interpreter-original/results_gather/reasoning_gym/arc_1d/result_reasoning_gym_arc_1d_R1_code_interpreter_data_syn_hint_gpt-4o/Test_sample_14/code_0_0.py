def transform_grid(input_grid):
    # Find the first and last non-zero indices
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    last_non_zero = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero:last_non_zero + 1]
    
    # Create the output grid
    output_grid = non_zero_sequence + [0] * (len(input_grid) - len(non_zero_sequence))
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
output = transform_grid(test_input)
print(output)