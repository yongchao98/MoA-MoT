def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def mark_rectangle(r, c):
        # Find the right and bottom boundaries of the rectangle
        right = c
        while right < cols and grid[r][right] in ('#', '█'):
            right += 1
        bottom = r
        while bottom < rows and grid[bottom][c] in ('#', '█'):
            bottom += 1
        # Mark the entire rectangle as visited
        for i in range(r, bottom):
            for j in range(c, right):
                visited[i][j] = True

    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in ('#', '█') and not visited[i][j]:
                mark_rectangle(i, j)
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                    ############                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    #          #                                ",
    "                                    ############                                "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)