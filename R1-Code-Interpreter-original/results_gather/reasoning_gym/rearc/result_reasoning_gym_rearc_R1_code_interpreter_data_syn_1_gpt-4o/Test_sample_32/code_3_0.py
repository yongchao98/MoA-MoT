def apply_new_hypothesis(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid
    # Check for a specific number and propagate it
    for i in range(len(output_grid)):
        for j in range(len(output_grid[i])):
            if output_grid[i][j] == 1 and j < len(output_grid[i]) - 1:
                output_grid[i][j + 1] = 1
    return output_grid

# Test input grid
input_grid = [
    [4, 1, 1, 0, 1, 1, 1],
    [1, 4, 1, 0, 1, 1, 1],
    [1, 1, 1, 0, 1, 1, 1]
]

# Apply the new hypothesis
output_grid = apply_new_hypothesis(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))