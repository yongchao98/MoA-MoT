def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    # Apply the transformation rule
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] in [0, 3, 5, 8, 9]:
                # Create a pattern around these key numbers
                if (i + j) % 2 == 0:
                    output_grid[i][j] = 0
                else:
                    output_grid[i][j] = 9
            else:
                # Default transformation for other numbers
                output_grid[i][j] = input_grid[i][j]

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