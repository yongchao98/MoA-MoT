# Define the input grid
input_grid = [
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 0, 1, 0, 0, 0, 2],
    [2, 2, 2, 2, 1, 0, 1, 2, 2, 2],
    [2, 2, 2, 2, 0, 1, 0, 0, 0, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
]

# Initialize the output grid with the same values as the input grid
output_grid = [row[:] for row in input_grid]

# Apply the mirroring rule
# Mirror 1s
output_grid[2][9] = 1
output_grid[6][5] = 1
output_grid[4][8] = 1
output_grid[8][4] = 1
output_grid[8][6] = 1

# Mirror 0s
output_grid[3][8] = 0
output_grid[7][4] = 0
output_grid[7][6] = 0
output_grid[7][7] = 0
output_grid[5][8] = 0
output_grid[9][4] = 0
output_grid[9][6] = 0
output_grid[9][7] = 0

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))