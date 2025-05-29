def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[8 for _ in range(cols)] for _ in range(rows)]

    # Check if the first row should be entirely 4
    if any(4 in row for row in input_grid):
        output_grid[0] = [4] * cols

    # Process each cell in the grid
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 2:
                # Set the current position to 2
                output_grid[r][c] = 2
                # Propagate 4 around the 2
                if r > 0:
                    output_grid[r-1][c] = 4
                if r < rows - 1:
                    output_grid[r+1][c] = 4
                if c > 0:
                    output_grid[r][c-1] = 4
                if c < cols - 1:
                    output_grid[r][c+1] = 4

    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 8],
    [8, 2, 8, 8, 8, 8, 8, 8, 8, 8, 2],
    [2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 2, 8, 8, 4, 8, 2, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 2, 8, 8, 8, 8, 2, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 2]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))