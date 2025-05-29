# Define the test input grid
input_grid = [
    [1, 9, 1, 1, 1, 1],
    [1, 1, 1, 1, 2, 1],
    [1, 2, 1, 1, 2, 1],
    [1, 1, 1, 1, 2, 1],
    [1, 1, 1, 1, 2, 1],
    [1, 0, 1, 1, 1, 1]
]

# Extract the subgrid from the bottom-right corner
output_grid = [row[-2:] for row in input_grid[-4:]]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))