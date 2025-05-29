def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 9:
                # Apply transformations around the cluster of 9s
                if j > 0 and input_grid[i][j-1] != 9:  # Left
                    output_grid[i][j-1] = 5
                if j < cols - 1 and input_grid[i][j+1] != 9:  # Right
                    output_grid[i][j+1] = 4
                if i > 0 and input_grid[i-1][j] != 9:  # Above
                    output_grid[i-1][j] = 1
                if i < rows - 1 and input_grid[i+1][j] != 9:  # Below
                    output_grid[i+1][j] = 8

    return output_grid

# Test input grid
input_grid = [
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 9, 9, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 9, 9, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 9, 9, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 9, 9, 3, 3],
    [3, 9, 9, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 9, 9, 3, 3, 3, 9, 9, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 9, 9, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))