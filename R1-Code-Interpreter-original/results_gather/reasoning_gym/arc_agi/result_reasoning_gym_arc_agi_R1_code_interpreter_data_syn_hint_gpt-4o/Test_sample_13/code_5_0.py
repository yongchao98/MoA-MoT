def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    def is_surrounded_by_4s(region):
        for r, c in region:
            directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
            for dr, dc in directions:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    if input_grid[nr][nc] != 4:
                        return False
                else:
                    return False
        return True

    def flood_fill(r, c, visited):
        stack = [(r, c)]
        region = []
        while stack:
            cr, cc = stack.pop()
            if (cr, cc) in visited:
                continue
            visited.add((cr, cc))
            region.append((cr, cc))
            directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
            for dr, dc in directions:
                nr, nc = cr + dr, cc + dc
                if 0 <= nr < rows and 0 <= nc < cols and input_grid[nr][nc] == 8 and (nr, nc) not in visited:
                    stack.append((nr, nc))
        return region

    visited = set()
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] == 8 and (r, c) not in visited:
                region = flood_fill(r, c, visited)
                if is_surrounded_by_4s(region):
                    for rr, cc in region:
                        output_grid[rr][cc] = 5

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