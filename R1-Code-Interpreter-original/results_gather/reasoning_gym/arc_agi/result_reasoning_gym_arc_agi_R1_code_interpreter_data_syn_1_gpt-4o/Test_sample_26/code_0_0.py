def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 1:
                # Check horizontal clusters
                if c > 0 and input_grid[r][c - 1] == 8:
                    # Start of a cluster
                    start = c
                    while c < cols and input_grid[r][c] == 1:
                        c += 1
                    end = c - 1
                    if end - start >= 1:
                        output_grid[r][start] = 0
                        output_grid[r][start + 1] = 0
                        output_grid[r][end] = 0
                        output_grid[r][end - 1] = 0

    for c in range(cols):
        for r in range(rows):
            if input_grid[r][c] == 1:
                # Check vertical clusters
                if r > 0 and input_grid[r - 1][c] == 8:
                    # Start of a cluster
                    start = r
                    while r < rows and input_grid[r][c] == 1:
                        r += 1
                    end = r - 1
                    if end - start >= 1:
                        output_grid[start][c] = 0
                        output_grid[start + 1][c] = 0
                        output_grid[end][c] = 0
                        output_grid[end - 1][c] = 0

    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 1, 1],
    [8, 8, 8, 8, 8, 1, 1],
    [8, 8, 1, 1, 8, 8, 8],
    [8, 8, 1, 1, 8, 8, 8],
    [8, 8, 8, 8, 1, 1, 8],
    [8, 8, 8, 8, 1, 1, 8],
    [8, 8, 8, 8, 8, 8, 8],
    [8, 1, 1, 8, 8, 8, 8],
    [8, 1, 1, 8, 8, 1, 1],
    [8, 8, 8, 8, 8, 1, 1],
    [8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 1, 1, 8],
    [8, 8, 8, 8, 1, 1, 8],
    [8, 1, 1, 8, 8, 8, 8],
    [8, 1, 1, 8, 8, 8, 8]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))