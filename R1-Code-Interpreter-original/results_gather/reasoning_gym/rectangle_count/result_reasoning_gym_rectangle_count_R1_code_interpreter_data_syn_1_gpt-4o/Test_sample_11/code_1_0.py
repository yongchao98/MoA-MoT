def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if grid[r][c] not in ('#', '█') or visited[r][c]:
            return False
        # Mark the current position as visited
        visited[r][c] = True
        # Check the boundaries of the rectangle
        right = c
        while right < cols and grid[r][right] in ('#', '█'):
            right += 1
        bottom = r
        while bottom < rows and grid[bottom][c] in ('#', '█'):
            bottom += 1
        # Check if the rectangle is valid
        for i in range(r, bottom):
            for j in range(c, right):
                if grid[i][j] not in ('#', '█'):
                    return False
                visited[i][j] = True
        return True

    for i in range(rows):
        for j in range(cols):
            if is_rectangle(i, j):
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