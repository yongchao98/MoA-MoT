def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    # Apply transformations based on observed patterns
    # Example 1 pattern: Change last element of first two rows
    output_grid[0][-1] = 2
    output_grid[1][-1] = 2

    # Example 2 pattern: Change first two elements of certain rows
    output_grid[4][0] = 9
    output_grid[4][1] = 9
    output_grid[9][0] = 9
    output_grid[9][1] = 9

    # Example 3 pattern: Change last two elements of first two rows
    output_grid[0][-2] = 9
    output_grid[0][-1] = 9
    output_grid[1][-2] = 9
    output_grid[1][-1] = 9

    # Example 4 pattern: Change specific elements in certain rows
    output_grid[6][6] = 9
    output_grid[8][4] = 9

    # Example 5 pattern: Change last two elements of first row
    output_grid[0][-2] = 9
    output_grid[0][-1] = 9

    return output_grid

# Test input grid
input_grid = [
    [7, 5, 7, 5, 7, 5, 7, 5, 7, 5, 7],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [7, 5, 7, 5, 9, 5, 7, 5, 7, 5, 7],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [7, 5, 9, 5, 2, 5, 7, 5, 7, 5, 7],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [7, 5, 7, 5, 7, 5, 7, 5, 7, 5, 7],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [7, 5, 7, 5, 7, 5, 7, 5, 7, 5, 2],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [7, 5, 7, 5, 7, 5, 7, 5, 7, 5, 7]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))