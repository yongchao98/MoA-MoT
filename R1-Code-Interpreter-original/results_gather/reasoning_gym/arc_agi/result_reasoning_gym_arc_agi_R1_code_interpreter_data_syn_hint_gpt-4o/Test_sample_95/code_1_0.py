def transform_grid(input_grid):
    output_grid = []
    
    # First three rows are the same
    for i in range(3):
        output_grid.append(input_grid[i])
    
    # Middle section transformation
    for i in range(3, 6):
        new_row = [min(input_grid[i][0], input_grid[i][1]), input_grid[i][1], input_grid[i][2]]
        if new_row[0] == 7:
            new_row[0] = input_grid[i][1]
        output_grid.append(new_row)
    
    # Last section transformation
    for i in range(6, 9):
        new_row = [min(input_grid[i][0], input_grid[i][1]), input_grid[i][1], input_grid[i][2]]
        if new_row[0] == 7:
            new_row[0] = input_grid[i][1]
        output_grid.append(new_row)
    
    return output_grid

# Test input grid
input_grid = [
    [7, 0, 0],
    [7, 7, 9],
    [7, 7, 7],
    [7, 5, 9],
    [7, 9, 7],
    [7, 7, 7],
    [9, 6, 6],
    [7, 7, 9],
    [7, 7, 7],
    [7, 4, 9],
    [7, 4, 7]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))