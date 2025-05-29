def transform_grid(input_grid):
    # Find the non-zero sequence
    non_zero_sequence = [x for x in input_grid if x != 0]
    
    # Create the output grid with zeros
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence shifted by 3 positions
    for i, value in enumerate(non_zero_sequence):
        output_grid[i + 3] = value
    
    return output_grid

# Test input
input_grid = [0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)