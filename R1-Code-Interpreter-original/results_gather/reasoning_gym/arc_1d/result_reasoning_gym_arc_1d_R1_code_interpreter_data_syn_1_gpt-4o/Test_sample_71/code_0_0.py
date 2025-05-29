def transform_grid(input_grid):
    # Find the position of the first non-zero element from the end
    non_zero_sequence = []
    for i in range(len(input_grid) - 1, -1, -1):
        if input_grid[i] != 0:
            non_zero_sequence.insert(0, input_grid[i])
        else:
            break
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence starting from the 6th position
    for i in range(len(non_zero_sequence)):
        output_grid[5 + i] = non_zero_sequence[i]
    
    # Find the position of '2' in the input grid and place it in the 8th position
    if 2 in input_grid:
        output_grid[7] = 2
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 2, 0, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)