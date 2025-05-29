def transform_grid(input_grid):
    # Create the output grid
    output_grid = [0] * len(input_grid)
    
    # Find the position of '2' in the input grid
    two_position = input_grid.index(2)
    
    # Place '2' in the 9th position of the output grid
    output_grid[8] = 2
    
    # Find the non-zero sequence after '2'
    non_zero_sequence = []
    for i in range(two_position + 1, len(input_grid)):
        if input_grid[i] != 0:
            non_zero_sequence.append(input_grid[i])
    
    # Place the non-zero sequence starting from the 6th position
    for i in range(len(non_zero_sequence)):
        output_grid[5 + i] = non_zero_sequence[i]
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 2, 0, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)