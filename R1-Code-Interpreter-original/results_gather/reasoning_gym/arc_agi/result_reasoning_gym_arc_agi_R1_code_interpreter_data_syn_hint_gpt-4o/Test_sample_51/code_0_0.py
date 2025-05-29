def repeat_grid(input_grid):
    # Get the dimensions of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create the output grid by repeating the input grid twice in both dimensions
    output_grid = []
    for i in range(rows):
        new_row = []
        for j in range(cols):
            new_row.extend([input_grid[i][j]] * 2)
        output_grid.append(new_row)
        output_grid.append(new_row)
    
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 3, 3],
    [9, 3, 9, 9],
    [9, 9, 9, 3],
    [3, 9, 3, 9]
]

# Generate the output grid
output_grid = repeat_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))