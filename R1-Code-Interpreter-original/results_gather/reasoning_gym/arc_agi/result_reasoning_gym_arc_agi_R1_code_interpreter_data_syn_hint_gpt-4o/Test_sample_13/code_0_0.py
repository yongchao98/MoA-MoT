def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    def is_surrounded_by_4s(r, c):
        # Check if the cell at (r, c) is surrounded by 4s or edges
        directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                if input_grid[nr][nc] != 4:
                    return False
        return True

    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 8 and is_surrounded_by_4s(r, c):
                output_grid[r][c] = 5

    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 4, 4, 4, 4, 4],
    [8, 8, 8, 8, 4, 4, 8, 4, 4, 8, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 8, 4, 4, 4, 4, 4, 8, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 8, 8, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 4, 4],
    [8, 8, 8, 8, 4, 4, 8, 8, 8, 4, 4, 4, 4, 8, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8],
    [4, 4, 4, 4, 4, 8, 8, 8, 8, 4, 4, 4, 4, 4, 4, 8]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))