# Define the input grid
input_grid = [
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 8, 3, 3, 3, 3, 3, 3, 1, 3],
    [3, 3, 9, 9, 9, 9, 9, 9, 3, 3],
    [3, 3, 9, 3, 3, 3, 3, 9, 3, 3],
    [3, 3, 9, 3, 3, 3, 3, 9, 3, 3],
    [3, 3, 9, 3, 3, 3, 3, 9, 3, 3],
    [3, 3, 9, 3, 3, 3, 3, 9, 3, 3],
    [3, 3, 9, 9, 9, 9, 9, 9, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
]

# Create a copy of the input grid for the output
output_grid = [row[:] for row in input_grid]

# Apply the transformation rule
output_grid[3][3] = output_grid[3][4] = output_grid[4][3] = output_grid[4][4] = input_grid[1][1]
output_grid[3][5] = output_grid[3][6] = output_grid[4][5] = output_grid[4][6] = input_grid[1][8]
output_grid[5][3] = output_grid[5][4] = output_grid[6][3] = output_grid[6][4] = input_grid[1][8]
output_grid[5][5] = output_grid[5][6] = output_grid[6][5] = output_grid[6][6] = input_grid[1][1]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))