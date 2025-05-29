def transform_grid(input_grid):
    # Find the non-zero sequence
    non_zero_sequence = [x for x in input_grid if x != 0]
    
    # Create the output grid with the non-zero sequence at the start
    output_grid = non_zero_sequence + [0] * (len(input_grid) - len(non_zero_sequence))
    
    return output_grid

# Test input grid
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
print(output_grid)