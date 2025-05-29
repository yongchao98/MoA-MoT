def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if grid[r][c] != '#' or visited[r][c]:
            return False
        # Find the bottom-right corner of the rectangle
        start_r, start_c = r, c
        end_r, end_c = r, c
        while end_c < cols and grid[start_r][end_c] == '#':
            end_c += 1
        end_c -= 1
        while end_r < rows and grid[end_r][start_c] == '#':
            end_r += 1
        end_r -= 1
        # Check if the rectangle is closed
        for i in range(start_r, end_r + 1):
            if grid[i][start_c] != '#' or grid[i][end_c] != '#':
                return False
        for j in range(start_c, end_c + 1):
            if grid[start_r][j] != '#' or grid[end_r][j] != '#':
                return False
        # Mark the rectangle as visited
        for i in range(start_r, end_r + 1):
            for j in range(start_c, end_c + 1):
                visited[i][j] = True
        return True

    for r in range(rows):
        for c in range(cols):
            if is_rectangle(r, c):
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                               #####################            ",
    "                                               #                   #            ",
    "                                               #####################            "
]

# Count the rectangles
print(count_rectangles(grid))