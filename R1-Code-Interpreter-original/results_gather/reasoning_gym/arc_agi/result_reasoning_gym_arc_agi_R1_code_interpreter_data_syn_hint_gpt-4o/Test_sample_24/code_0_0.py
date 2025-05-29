def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Start with a copy of the input grid

    # Function to propagate '3's horizontally and vertically
    def propagate_3s(r, c):
        # Propagate horizontally
        for i in range(c, cols):
            if output_grid[r][i] == 7:
                output_grid[r][i] = 3
            else:
                break
        for i in range(c, -1, -1):
            if output_grid[r][i] == 7:
                output_grid[r][i] = 3
            else:
                break
        # Propagate vertically
        for i in range(r, rows):
            if output_grid[i][c] == 7:
                output_grid[i][c] = 3
            else:
                break
        for i in range(r, -1, -1):
            if output_grid[i][c] == 7:
                output_grid[i][c] = 3
            else:
                break

    # Identify starting points for '3's
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 3:
                propagate_3s(r, c)
            elif input_grid[r][c] == 1:
                # Check if this '1' should start a propagation
                if (r > 0 and input_grid[r-1][c] == 3) or (c > 0 and input_grid[r][c-1] == 3):
                    propagate_3s(r, c)

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
for row in output_grid:
    print(' '.join(map(str, row)))