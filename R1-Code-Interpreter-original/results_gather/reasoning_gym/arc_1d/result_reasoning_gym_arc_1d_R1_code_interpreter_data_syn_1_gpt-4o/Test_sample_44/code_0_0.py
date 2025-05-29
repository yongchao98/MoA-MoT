def transform_grid(input_grid):
    # Find the position of the number 2
    pos_2 = input_grid.index(2)
    
    # Find the first non-zero number after the 2
    for i in range(pos_2 + 1, len(input_grid)):
        if input_grid[i] != 0:
            start_pos = i
            break
    
    # Create the output grid
    output_grid = input_grid[:pos_2] + input_grid[start_pos:] + [0] * (start_pos - pos_2)
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)