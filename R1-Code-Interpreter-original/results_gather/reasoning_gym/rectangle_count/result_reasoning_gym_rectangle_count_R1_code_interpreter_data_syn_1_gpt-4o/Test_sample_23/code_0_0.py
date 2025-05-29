def count_rectangles(grid):
    def is_top_left_corner(x, y):
        return (grid[x][y] == '#' and
                (x == 0 or grid[x-1][y] == ' ') and
                (y == 0 or grid[x][y-1] == ' '))

    def mark_rectangle(x, y):
        # Mark the rectangle as visited
        max_x, max_y = x, y
        while max_x < len(grid) and grid[max_x][y] == '#':
            max_x += 1
        while max_y < len(grid[0]) and grid[x][max_y] == '#':
            max_y += 1
        for i in range(x, max_x):
            for j in range(y, max_y):
                visited[i][j] = True

    visited = [[False] * len(grid[0]) for _ in range(len(grid))]
    rectangle_count = 0

    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if is_top_left_corner(i, j) and not visited[i][j]:
                rectangle_count += 1
                mark_rectangle(i, j)

    return rectangle_count

# Define the grid
grid = [
    "                       ######################                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       # ############       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # #          #       #                                   ",
    "                       # ############       #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       #                    #                                   ",
    "                       ######################                                   "
]

# Count the rectangles
print(count_rectangles(grid))