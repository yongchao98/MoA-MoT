def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        if grid[r][c] != '#' and grid[r][c] != '█':
            return False
        if visited[r][c]:
            return False
        return True

    def mark_rectangle(r, c):
        nonlocal rectangle_count
        # Find the bottom-right corner
        end_r, end_c = r, c
        while end_c < cols and (grid[r][end_c] == '#' or grid[r][end_c] == '█'):
            end_c += 1
        end_c -= 1
        while end_r < rows and (grid[end_r][c] == '#' or grid[end_r][c] == '█'):
            end_r += 1
        end_r -= 1

        # Mark the rectangle as visited
        for i in range(r, end_r + 1):
            for j in range(c, end_c + 1):
                visited[i][j] = True

        rectangle_count += 1

    for i in range(rows):
        for j in range(cols):
            if is_rectangle(i, j):
                mark_rectangle(i, j)

    return rectangle_count

# Define the grid
grid = [
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "  #########################################                                      ",
    "  #                                       #                                      ",
    "  #                                       #                                      ",
    "  #########################################                                      ",
    "                                       #########                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                                 ",
    "                                       #       #                             ### ",
    "                                       #########                             # # ",
    "                                                                             # # ",
    "                                                                             # # ",
    "                                                                             # # ",
    "                                                                             # # ",
    "                                                                             # # ",
    "                                                                             # # ",
    "                                                                             # # ",
    "                                                                             # # ",
    "                                                                             # # ",
    "                                                                             ### ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)