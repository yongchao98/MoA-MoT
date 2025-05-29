def transform_grid(input_grid):
    # Find the first non-zero sequence from the left
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    first_non_zero_value = input_grid[first_non_zero_index]
    
    # Find the last non-zero sequence from the right
    last_non_zero_index = next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)
    last_non_zero_value = input_grid[-(last_non_zero_index + 1)]
    
    # Create the output grid
    output_grid = input_grid[:]
    
    # Replace leading zeros with the first non-zero value
    for i in range(first_non_zero_index):
        output_grid[i] = first_non_zero_value
    
    # Replace trailing zeros with the last non-zero value
    for i in range(len(input_grid) - last_non_zero_index, len(input_grid)):
        output_grid[i] = last_non_zero_value
    
    return output_grid

# Test input
input_grid = [0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)