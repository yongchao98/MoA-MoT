# Define the input grid
input_grid = [
    [9, 6, 3, 6, 3, 6],
    [6, 6, 6, 6, 6, 6],
    [9, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6],
    [1, 6, 0, 6, 0, 6],
    [6, 6, 6, 6, 6, 6]
]

# Extract the 3x3 grid
output_grid = [
    [input_grid[i][j] for j in [0, 2, 4]]
    for i in [0, 2, 4]
]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))