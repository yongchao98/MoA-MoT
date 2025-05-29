def transform_grid(input_grid):
    # Find the first non-zero index
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero_index:]
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    output_grid[first_non_zero_index - 3:first_non_zero_index - 3 + len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)