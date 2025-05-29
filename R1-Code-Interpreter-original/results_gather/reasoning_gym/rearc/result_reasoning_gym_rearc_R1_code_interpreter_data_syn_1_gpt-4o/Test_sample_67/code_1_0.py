def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy input grid to output grid

    # Find the block of zeros to transform
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if input_grid[i][j] == 0:
                # Check if this zero is part of a block surrounded by non-zero numbers
                if (input_grid[i-1][j] != 0 and input_grid[i+1][j] != 0 and
                    input_grid[i][j-1] != 0 and input_grid[i][j+1] != 0):
                    output_grid[i][j] = 1

    return output_grid

# Test input grid
input_grid = [
    [3, 3, 0, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 0],
    [3, 3, 3, 3, 0, 0, 0],
    [3, 0, 0, 3, 0, 0, 0],
    [0, 3, 0, 3, 0, 0, 0],
    [0, 3, 0, 0, 0, 3, 3],
    [3, 3, 0, 3, 0, 3, 3]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))