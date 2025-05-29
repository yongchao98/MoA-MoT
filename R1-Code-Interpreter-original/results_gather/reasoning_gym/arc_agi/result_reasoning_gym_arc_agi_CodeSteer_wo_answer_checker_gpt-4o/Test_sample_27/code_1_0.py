def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    rows = len(input_grid)
    cols = len(input_grid[0])

    for i in range(rows - 1):
        for j in range(cols - 1):
            # Check for the top-left corner of a '9' block
            if input_grid[i][j] == 9 and input_grid[i][j + 1] == 9 and input_grid[i + 1][j] == 9 and input_grid[i + 1][j + 1] == 9:
                # Place '5' to the left of the top-left '9'
                if j > 0 and output_grid[i][j - 1] == 3:
                    output_grid[i][j - 1] = 5
                # Place '4' to the right of the top-right '9'
                if j + 2 < cols and output_grid[i][j + 2] == 3:
                    output_grid[i][j + 2] = 4
                # Place '1' below the bottom-left '9'
                if i + 2 < rows and output_grid[i + 2][j] == 3:
                    output_grid[i + 2][j] = 1
                # Place '8' below the bottom-right '9'
                if i + 2 < rows and j + 1 < cols and output_grid[i + 2][j + 1] == 3:
                    output_grid[i + 2][j + 1] = 8

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