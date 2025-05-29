# Test input grid
input_grid = [
    [6, 4, 6, 4, 6, 0],
    [6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 0],
    [6, 6, 6, 6, 6, 6],
    [6, 1, 6, 1, 6, 7],
    [6, 6, 6, 6, 6, 6]
]

# Extracting the output grid based on the identified pattern
output_grid = [
    [input_grid[0][1], input_grid[0][3], input_grid[0][5]],
    [input_grid[2][1], input_grid[2][3], input_grid[2][5]],
    [input_grid[4][1], input_grid[4][3], input_grid[4][5]]
]

# Printing the output grid
for row in output_grid:
    print(' '.join(map(str, row)))