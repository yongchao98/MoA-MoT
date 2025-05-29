def max_area_of_island(grid):
    def dfs(x, y):
        if x < 0 or x >= len(grid) or y < 0 or y >= len(grid[0]) or grid[x][y] == 0 or visited[x][y]:
            return 0
        visited[x][y] = True
        area = 1  # Current cell
        # Explore neighbors
        area += dfs(x + 1, y)
        area += dfs(x - 1, y)
        area += dfs(x, y + 1)
        area += dfs(x, y - 1)
        return area

    if not grid or not grid[0]:
        return 0

    max_area = 0
    visited = [[False for _ in range(len(grid[0]))] for _ in range(len(grid))]

    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == 1 and not visited[i][j]:
                max_area = max(max_area, dfs(i, j))

    return max_area

# Define the grid
grid = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
]

# Calculate the maximum area of an island
result = max_area_of_island(grid)
print(result)