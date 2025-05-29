def count_rectangles(grid):
    def is_corner(x, y):
        return (grid[x][y] == '#' and
                ((x > 0 and grid[x-1][y] == ' ') or x == 0) and
                ((y > 0 and grid[x][y-1] == ' ') or y == 0))

    def find_rectangle(x, y):
        max_x, max_y = x, y
        while max_x + 1 < len(grid) and grid[max_x + 1][y] in ['#', '█']:
            max_x += 1
        while max_y + 1 < len(grid[0]) and grid[x][max_y + 1] in ['#', '█']:
            max_y += 1
        return max_x, max_y

    visited = set()
    rectangle_count = 0

    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if is_corner(i, j) and (i, j) not in visited:
                max_x, max_y = find_rectangle(i, j)
                for x in range(i, max_x + 1):
                    for y in range(j, max_y + 1):
                        visited.add((x, y))
                rectangle_count += 1

    return rectangle_count

grid = [
    "                                      ########################                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      #                      #                  ",
    "                                      ########################                  ",
    "                                                             ###                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                     ########█#█#######         ",
    "                                                     #       # #      #         ",
    "        ####################                         #       # #      #         ",
    "        #                  #                         #       # #      #         ",
    "        #                  #                         #       # #      #         ",
    "        #                  #                         #       # #      #         ",
    "        #                  #                         #       # #      #         ",
    "        #                  #                         ########█#█#######         ",
    "        #                  #                                 # #                ",
    "        #                  #                                 # #                ",
    "        #                  #                                 # #                ",
    "        #                  #                                 # #                ",
    "        #                  #                                 # #                ",
    "        ####################                                 # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             ###                ",
    "                           ###################################################  ",
    "                           #                                                 #  ",
    "                           #                                                 #  ",
    "                           #                                                 #  ",
    "                           ###################################################  ",
]

print(count_rectangles(grid))