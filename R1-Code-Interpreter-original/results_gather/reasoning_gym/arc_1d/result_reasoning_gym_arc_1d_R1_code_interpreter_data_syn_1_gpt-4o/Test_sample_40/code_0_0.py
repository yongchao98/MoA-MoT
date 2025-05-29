def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if first_non_zero is None:
        return input_grid  # No non-zero elements, return as is
    
    # Find the last non-zero element
    last_non_zero = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero:last_non_zero + 1]
    
    # Create the new grid with the sequence shifted two positions to the left
    new_grid = [0] * len(input_grid)
    new_grid[2:2 + len(non_zero_sequence)] = non_zero_sequence
    
    return new_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)