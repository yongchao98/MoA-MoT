def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = [[False] * cols for _ in range(rows)]
    rectangle_count = 0

    def is_rectangle(r, c):
        return grid[r][c] == '#' or grid[r][c] == '█'

    def mark_rectangle(r, c):
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_c < cols and is_rectangle(r, end_c):
            end_c += 1
        while end_r < rows and is_rectangle(end_r, c):
            end_r += 1
        end_r -= 1
        end_c -= 1

        # Mark the rectangle as visited
        for i in range(r, end_r + 1):
            for j in range(c, end_c + 1):
                visited[i][j] = True

    for r in range(rows):
        for c in range(cols):
            if is_rectangle(r, c) and not visited[r][c]:
                mark_rectangle(r, c)
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                       #############                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                       #           #                                            ",
    "                    ###█###########█######################                      ",
    "                    #  #############                     #                      ",
    "                    #                                    #                      ",
    "                    #                                    #                      ",
    "                    #                                    #                      ",
    "                    #                                    #                      ",
    "                    #                     ###############█#########             ",
    "                    #                     #              #        #             ",
    "                    #                     #              #        #             ",
    "                    #                     ###############█#########             ",
    "                    #                                    #                      ",
    "                    #                                    #                      ",
    "                    #                                    #                      ",
    "                    #                                    #                      ",
    "                    #                                    #                      ",
    "               #####█####################################█#####                 ",
    "               #    #                                    #    #                 ",
    "               #    #                                    #    #                 ",
    "               #    #                                    #    #                 ",
    "               #    #                                    #    #                 ",
    "               #    #                                    #    #                 ",
    "               #    #                                    #    #                 ",
    "               #    #                                    #    #                 ",
    "               #    #                                    #    #                 ",
    "               #    ######################################    #                 ",
    "               #                                              #                 ",
    "               ################################################                 ",
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)