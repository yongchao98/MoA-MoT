# Define the input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 7, 7],
    [5, 5, 5, 7, 7, 7, 7, 7, 8, 8, 7, 7],
    [7, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [5, 5, 5, 7, 7, 7, 3, 3, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 0, 7, 7, 7, 7]
]

# Extract the subgrid from the bottom-right corner
output_grid = [row[7:] for row in input_grid[8:]]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))