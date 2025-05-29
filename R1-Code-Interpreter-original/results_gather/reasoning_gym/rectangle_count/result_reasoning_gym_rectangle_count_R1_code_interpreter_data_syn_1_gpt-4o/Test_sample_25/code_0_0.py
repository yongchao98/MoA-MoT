def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if grid[r][c] != '#' or visited[r][c]:
            return False
        # Mark the entire rectangle as visited
        start_r, start_c = r, c
        while c < cols and grid[r][c] == '#':
            c += 1
        end_c = c - 1
        while r < rows and grid[r][start_c] == '#':
            r += 1
        end_r = r - 1
        for i in range(start_r, end_r + 1):
            for j in range(start_c, end_c + 1):
                if grid[i][j] != '#':
                    return False
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