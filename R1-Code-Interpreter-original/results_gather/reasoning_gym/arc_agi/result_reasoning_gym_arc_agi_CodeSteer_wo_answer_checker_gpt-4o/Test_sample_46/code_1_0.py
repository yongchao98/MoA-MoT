def flood_fill(grid, x, y, target, visited):
    rows, cols = len(grid), len(grid[0])
    stack = [(x, y)]
    component = []

    while stack:
        cx, cy = stack.pop()
        if (0 <= cx < rows and 0 <= cy < cols and
            grid[cx][cy] == target and not visited[cx][cy]):
            visited[cx][cy] = True
            component.append((cx, cy))
            # Add neighbors to stack
            stack.extend([(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)])
    
    return component

def is_boundary(cell, grid):
    x, y = cell
    rows, cols = len(grid), len(grid[0])
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if nx < 0 or nx >= rows or ny < 0 or ny >= cols or grid[nx][ny] == 7:
            return True
    return False

def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid
    visited = [[False] * cols for _ in range(rows)]

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != 7 and not visited[i][j]:
                # Find the component
                component = flood_fill(input_grid, i, j, input_grid[i][j], visited)
                # Determine the boundary of the component
                for cell in component:
                    if is_boundary(cell, input_grid):
                        x, y = cell
                        output_grid[x][y] = 7

    return output_grid

# Test input grid
input_grid = [
    [6, 7, 7, 7, 6, 7, 7, 7, 7, 6],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 6, 6, 6, 7, 6, 7, 6, 7],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [6, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 6, 6, 7, 7],
    [7, 6, 7, 7, 7, 7, 6, 6, 7, 7],
    [7, 7, 7, 7, 7, 7, 6, 6, 6, 7],
    [7, 7, 6, 7, 7, 6, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))