def count_rectangles(grid):
    def is_corner(x, y):
        # Check if the current position is a top-left corner of a rectangle
        if grid[x][y] in ['#', '█']:
            if (x == 0 or grid[x-1][y] == ' ') and (y == 0 or grid[x][y-1] == ' '):
                return True
        return False

    def find_rectangle(x, y):
        # Find the bottom-right corner of the rectangle starting from (x, y)
        max_x, max_y = x, y
        while max_x + 1 < len(grid) and grid[max_x + 1][y] in ['#', '█']:
            max_x += 1
        while max_y + 1 < len(grid[0]) and grid[x][max_y + 1] in ['#', '█']:
            max_y += 1
        return max_x, max_y

    visited = set()
    rectangle_count = 0

    for i in range(len(grid)):
        for j in range(len(grid[i])):  # Use len(grid[i]) to handle varying row lengths
            if is_corner(i, j) and (i, j) not in visited:
                bottom_right = find_rectangle(i, j)
                # Mark all positions of this rectangle as visited
                for x in range(i, bottom_right[0] + 1):
                    for y in range(j, bottom_right[1] + 1):
                        visited.add((x, y))
                rectangle_count += 1

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
    "       ###############                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #             #                                                           ",
    "       #  ###########█#######██████████████##################                   ",
    "       #  #          #       #            #                 #                   ",
    "       #  #          #       #            #                 #                   ",
    "       #  #          #       #            #                 #                   ",
    "       #  #          #       #            #                 #                   ",
    "       #  #          #       #            #                 #                   ",
    "       ###█###########       #            #                 #                   ",
    "          #                  #            #                 #                   ",
    "          #                  #            #                 #                   ",
    "          #                  #         ###█#################█################   ",
    "          #                  #         #  #                 #               #   ",
    "          #                  #         #  #                 #               #   ",
    "          #                  #         #  #                 #               #   ",
    "          #                  #         #  #                 #               #   ",
    "          #                  #         #  #                 #               #   ",
    "          #                  #         #  #                 #  #######      #   ",
    "          #                  #         ###█#################█##█#####█#######   ",
    "          #                  #            #                 #  #     #          ",
    "          ###################█#############                 #  #     #          ",
    "                             ################################  #     #          ",
    "                                                               #     #          ",
    "                                                               #######          ",
    "                                                                                 "
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)