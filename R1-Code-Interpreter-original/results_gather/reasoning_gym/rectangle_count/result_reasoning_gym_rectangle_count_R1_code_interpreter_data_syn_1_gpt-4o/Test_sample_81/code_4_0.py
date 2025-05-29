def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if grid[r][c] != '#' or visited[r][c]:
            return False
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_c < cols and grid[r][end_c] == '#':
            end_c += 1
        while end_r < rows and grid[end_r][c] == '#':
            end_r += 1
        # Check if the found area is a rectangle
        for i in range(r, end_r):
            for j in range(c, end_c):
                if grid[i][j] != '#' or visited[i][j]:
                    return False
                visited[i][j] = True
        # Ensure the rectangle is closed
        for i in range(r, end_r):
            if grid[i][end_c - 1] != '#':
                return False
        for j in range(c, end_c):
            if grid[end_r - 1][j] != '#':
                return False
        return True

    for i in range(rows):
        for j in range(cols):
            if is_rectangle(i, j):
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                   ####                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   #  #                                         ",
    "                                   ####                                         "
]

# Count rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)