def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 1:
                # Check all surrounding positions and replace '7' with '3'
                for di in [-1, 0, 1]:
                    for dj in [-1, 0, 1]:
                        ni, nj = i + di, j + dj
                        while 0 <= ni < rows and 0 <= nj < cols and input_grid[ni][nj] == 7:
                            output_grid[ni][nj] = 3
                            ni += di
                            nj += dj

    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 1, 7, 7, 7, 1, 7, 7, 7, 1, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 1, 1, 1, 7, 7],
    [7, 7, 7, 1, 1, 7, 1, 7, 7, 7, 7, 1, 7],
    [1, 7, 7, 7, 1, 7, 1, 7, 7, 1, 7, 7, 3],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 1, 7],
    [1, 1, 7, 7, 1, 7, 1, 1, 7, 7, 7, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 1, 7, 7, 7, 1, 7],
    [7, 1, 7, 1, 7, 7, 1, 7, 7, 1, 7, 7, 3],
    [7, 7, 7, 7, 7, 1, 7, 1, 7, 7, 7, 7, 7],
    [7, 1, 7, 7, 7, 7, 7, 7, 1, 1, 1, 7, 3],
    [7, 1, 7, 7, 7, 1, 7, 7, 1, 7, 1, 7, 1],
    [1, 1, 7, 7, 7, 7, 7, 7, 7, 1, 7, 7, 1],
    [1, 7, 7, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))