def count_rectangles(grid):
    def is_top_left_corner(x, y):
        # Check if the current position is a top-left corner of a rectangle
        if grid[x][y] not in ['#', '█']:
            return False
        # Check if it's the start of a horizontal and vertical line
        if (x == 0 or grid[x-1][y] == ' ') and (y == 0 or grid[x][y-1] == ' '):
            return True
        return False

    def find_rectangle(x, y):
        # Find the bottom-right corner of the rectangle starting from (x, y)
        max_x, max_y = x, y
        while max_x < len(grid) and grid[max_x][y] in ['#', '█']:
            max_x += 1
        while max_y < len(grid[0]) and grid[x][max_y] in ['#', '█']:
            max_y += 1
        return max_x - 1, max_y - 1

    visited = set()
    rectangle_count = 0

    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if (i, j) not in visited and is_top_left_corner(i, j):
                bottom_right_x, bottom_right_y = find_rectangle(i, j)
                # Mark all positions of this rectangle as visited
                for x in range(i, bottom_right_x + 1):
                    for y in range(j, bottom_right_y + 1):
                        visited.add((x, y))
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                          ##################################### ",
    "                                          #                                   # ",
    "                                          ##################################### ",
    "                                                           ###############      ",
    "                                                           #             #      ",
    "                                                           #             #      ",
    "                                                           #             #      ",
    "                                                           #             #      ",
    "      #####################################################█##           #      ",
    "      #                                                    # #           #      ",
    "      #                                                    # #           #      ",
    "      #                                                    # #           #      ",
    "      #                                                    # #           #      ",
    "      #                                                    # #           #      ",
    "      #                                                    # #           #      ",
    "      #                                                    # #           #      ",
    "      #   #############################################    # #           #      ",
    "      #   #                                           #    # #           #      ",
    "      #   #                                           #    # #           #      ",
    "      #   #                                           #    # #           #      ",
    "      #   #                                           #    # #           #      ",
    "      #   #                                           #    # #           #      ",
    "      #   #                                           #    # #           #      ",
    "      #   #                                           #    # #           #      ",
    "      #   #                                           #    ##█############      ",
    "      #   #                                           #      #                  ",
    "      #   #                                           #      #                  ",
    "      #   #                                           #      #                  ",
    "      ####█###########################################█#######                  ",
    "          #                                           #                         ",
    "          #                                           #                         ",
    "          #                                           #                         ",
    "          #############################################                         ",
    "                                      ###############################           ",
    "                                      #                             #           ",
    "                                      #                             #           ",
    "                                      #                             #           ",
    "                                      #                             #           ",
    "                                      #                             #           ",
    "                                      #                             #           ",
    "                                      #                             #           ",
    "                                      #                             #           ",
    "                                      ###############################           ",
    "                                                    #####################       ",
    "                                                    #                   #       ",
    "                                                    #####################       "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)