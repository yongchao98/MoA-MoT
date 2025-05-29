from collections import deque

def find_blocks(grid):
    rows, cols = len(grid), len(grid[0])
    visited = [[False] * cols for _ in range(rows)]
    blocks = []

    def bfs(x, y, value):
        queue = deque([(x, y)])
        block = []
        while queue:
            cx, cy = queue.popleft()
            if (0 <= cx < rows and 0 <= cy < cols and
                grid[cx][cy] == value and not visited[cx][cy]):
                visited[cx][cy] = True
                block.append((cx, cy))
                # Add neighbors to queue
                queue.extend([(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)])
        return block

    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 7 and not visited[i][j]:
                block = bfs(i, j, grid[i][j])
                blocks.append(block)

    return blocks

def is_boundary(cell, grid):
    x, y = cell
    rows, cols = len(grid), len(grid[0])
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if nx < 0 or nx >= rows or ny < 0 or ny >= cols or grid[nx][ny] == 7:
            return True
    return False

def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid
    blocks = find_blocks(input_grid)

    for block in blocks:
        for cell in block:
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