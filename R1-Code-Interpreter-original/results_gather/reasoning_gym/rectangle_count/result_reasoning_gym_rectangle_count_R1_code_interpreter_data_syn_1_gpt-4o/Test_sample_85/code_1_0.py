def count_rectangles(grid):
    # Convert each row of the grid to a list of characters
    grid = [list(row) for row in grid]
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = set()

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in ['#', '█']:
                # Find the bottom-right corner of the rectangle
                br_r, br_c = r, c
                while br_r < rows and grid[br_r][c] in ['#', '█']:
                    br_r += 1
                while br_c < cols and grid[r][br_c] in ['#', '█']:
                    br_c += 1
                # Add the rectangle defined by top-left and bottom-right corners
                rectangles.add((r, c, br_r - 1, br_c - 1))
                # Mark the rectangle as visited
                for i in range(r, br_r):
                    for j in range(c, br_c):
                        grid[i][j] = '.'

    return len(rectangles)

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
total_rectangles = count_rectangles(grid)
print(total_rectangles)