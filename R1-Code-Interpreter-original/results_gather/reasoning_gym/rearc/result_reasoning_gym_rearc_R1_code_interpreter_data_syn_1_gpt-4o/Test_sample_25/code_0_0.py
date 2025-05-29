# Define the input grid
input_grid = [
    [7, 7, 8, 8, 8, 8, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7]
]

# Create the output grid by applying the rule
output_grid = [row[:] for row in input_grid]  # Copy the input grid

# Apply the transformation rule
for i in range(1, len(output_grid)):
    for j in range(2, 6):
        output_grid[i][j] = 7

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))