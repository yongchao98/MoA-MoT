def is_boundary(x, y, grid):
    """Check if the cell (x, y) is a boundary cell of '1's."""
    rows, cols = len(grid), len(grid[0])
    if grid[x][y] != 1:
        return False
    # Check if any neighboring cell is not '1'
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < rows and 0 <= ny < cols and grid[nx][ny] != 1:
            return True
    return False

def flood_fill(x, y, grid, visited):
    """Perform flood fill to find all connected '1's starting from (x, y)."""
    stack = [(x, y)]
    component = []
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        component.append((cx, cy))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] == 1 and (nx, ny) not in visited:
                stack.append((nx, ny))
    return component

def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid
    visited = set()

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 1 and (i, j) not in visited:
                component = flood_fill(i, j, input_grid, visited)
                # Transform interior '1's to '5's
                for x, y in component:
                    if not is_boundary(x, y, input_grid):
                        output_grid[x][y] = 5

    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1],
    [6, 6, 6, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1, 1, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6],
    [6, 6, 6, 1, 6, 6, 6, 1, 6, 6, 1, 6, 6, 6, 6],
    [1, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1, 6, 6, 6, 6],
    [1, 1, 1, 1, 6, 6, 6, 1, 6, 6, 6, 6, 6, 6, 6],
    [1, 6, 6, 1, 6, 6, 6, 1, 6, 6, 6, 1, 1, 1, 1],
    [1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 6, 6, 1],
    [1, 1, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 6, 6, 1],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 1, 1, 6, 6, 6, 1, 1, 1, 1, 1, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 6, 1, 6, 6, 6, 6, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 6, 1, 6, 6, 6, 6, 1, 6, 6, 6],
    [6, 1, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 6, 6, 6]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))