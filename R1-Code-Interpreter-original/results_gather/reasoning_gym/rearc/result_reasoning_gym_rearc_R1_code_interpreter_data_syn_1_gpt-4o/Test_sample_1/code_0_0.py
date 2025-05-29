def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    # Apply the transformation rule
    for i in range(rows):
        for j in range(cols):
            if (i + j) % 3 == 0:  # Example pattern: diagonal based on sum of indices
                output_grid[i][j] = 0  # Replace with 0 as observed in examples
            elif (i + j) % 3 == 1:
                output_grid[i][j] = 9  # Replace with 9 as observed in examples
            # Add more conditions if needed based on further analysis

    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9],
    [8, 2, 8, 8, 8, 8, 8, 8, 8, 9, 9],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 5, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))