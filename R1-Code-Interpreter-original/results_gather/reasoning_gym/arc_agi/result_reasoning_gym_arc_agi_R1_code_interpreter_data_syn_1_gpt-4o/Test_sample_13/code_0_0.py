def find_clusters_and_transform(grid):
    rows = len(grid)
    cols = len(grid[0])
    visited = [[False] * cols for _ in range(rows)]
    output_grid = [row[:] for row in grid]

    def dfs(r, c):
        stack = [(r, c)]
        cluster = []
        while stack:
            x, y = stack.pop()
            if 0 <= x < rows and 0 <= y < cols and not visited[x][y] and grid[x][y] == 8:
                visited[x][y] = True
                cluster.append((x, y))
                # Check all 8 directions
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]:
                    stack.append((x + dx, y + dy))
        return cluster

    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 8 and not visited[i][j]:
                cluster = dfs(i, j)
                if len(cluster) > 1:  # Only transform clusters larger than 1
                    for x, y in cluster:
                        output_grid[x][y] = 5

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

output_grid = find_clusters_and_transform(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))