def apply_correct_rule(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid
    # Propagate the number from the last column to the left in the first row
    output_grid[0][-2] = output_grid[0][-1]
    return output_grid

# Test input grid
input_grid = [
    [4, 1, 1, 0, 1, 1, 1],
    [1, 4, 1, 0, 1, 1, 1],
    [1, 1, 1, 0, 1, 1, 1]
]

# Apply the corrected rule
output_grid = apply_correct_rule(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))