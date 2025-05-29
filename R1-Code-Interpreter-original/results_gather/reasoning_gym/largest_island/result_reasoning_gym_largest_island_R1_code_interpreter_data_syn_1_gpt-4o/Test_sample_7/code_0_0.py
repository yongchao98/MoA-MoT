def max_area_of_island(grid):
    def dfs(x, y):
        if x < 0 or x >= len(grid) or y < 0 or y >= len(grid[0]) or grid[x][y] == 0:
            return 0
        grid[x][y] = 0  # Mark the cell as visited
        area = 1  # Current cell
        # Explore all 4 directions
        area += dfs(x + 1, y)
        area += dfs(x - 1, y)
        area += dfs(x, y + 1)
        area += dfs(x, y - 1)
        return area

    max_area = 0
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == 1:
                max_area = max(max_area, dfs(i, j))
    return max_area

# Given 6x7 binary matrix grid
grid = [
    [0, 1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0],
    [0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 1]
]

# Calculate and print the maximum area of an island
print(max_area_of_island(grid))