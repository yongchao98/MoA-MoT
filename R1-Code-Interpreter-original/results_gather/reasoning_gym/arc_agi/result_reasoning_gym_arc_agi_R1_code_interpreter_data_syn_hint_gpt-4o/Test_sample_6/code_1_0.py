def transform_grid(input_grid):
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Determine the stopping point for replacement
    replace_until = 0
    for i in range(num_rows):
        if any(cell == 6 or cell == 9 for cell in input_grid[i]):
            replace_until = i
            break
    
    # Create the output grid
    output_grid = []
    for i in range(num_rows):
        if i < replace_until:
            output_grid.append([8] * num_cols)
        else:
            output_grid.append(input_grid[i])
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6, 8],
    [8, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 8],
    [8, 6, 9, 6, 6, 6, 6, 8, 8, 8, 8, 8],
    [8, 6, 6, 6, 6, 9, 6, 8, 6, 6, 6, 6],
    [8, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 9, 6],
    [8, 8, 8, 8, 6, 6, 6, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 6, 6, 9, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 6, 6, 6, 8, 8, 8, 8, 8],
    [6, 6, 6, 8, 9, 6, 6, 8, 8, 8, 8, 8],
    [9, 6, 6, 8, 8, 8, 8, 8, 6, 9, 6, 9],
    [6, 6, 6, 8, 8, 8, 8, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 9, 6, 6]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))