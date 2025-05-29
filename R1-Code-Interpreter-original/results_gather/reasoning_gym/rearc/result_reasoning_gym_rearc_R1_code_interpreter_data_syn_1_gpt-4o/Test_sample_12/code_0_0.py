# Define the input grid
input_grid = [
    [8, 8, 8, 0, 0, 0, 0],
    [8, 8, 8, 0, 0, 0, 8],
    [0, 8, 8, 0, 0, 0, 0],
    [8, 8, 8, 0, 8, 8, 8],
    [0, 8, 8, 8, 8, 8, 8],
    [0, 8, 8, 8, 0, 8, 8]
]

# Create a copy of the input grid for the output
output_grid = [row[:] for row in input_grid]

# Apply the rule: change a specific 3x3 subgrid to 1s
# Based on the examples, let's assume the subgrid is in the middle of the grid
for i in range(2, 5):
    for j in range(2, 5):
        output_grid[i][j] = 1

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))