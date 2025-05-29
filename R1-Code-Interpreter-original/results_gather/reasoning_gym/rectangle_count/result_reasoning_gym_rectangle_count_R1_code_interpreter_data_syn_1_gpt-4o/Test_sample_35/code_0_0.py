def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def mark_rectangle(r, c):
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_r < rows and grid[end_r][c] in ['#', '█']:
            end_r += 1
        while end_c < cols and grid[r][end_c] in ['#', '█']:
            end_c += 1
        end_r -= 1
        end_c -= 1

        # Mark the rectangle as visited
        for i in range(r, end_r + 1):
            for j in range(c, end_c + 1):
                visited[i][j] = True

    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in ['#', '█'] and not visited[i][j]:
                rectangle_count += 1
                mark_rectangle(i, j)

    return rectangle_count

# Define the grid
grid = [
    "                                  #########################################     ",
    "                                  #                                       #     ",
    "                                  #                                       #     ",
    "                                  #                                       #     ",
    "                                  #                                       #     ",
    "              ####################█##################################     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   #                                 #     #     ",
    "              #                   ##################################█######     ",
    "              #                                                     #           ",
    "              #                                                     #           ",
    "              #######################################################           "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)