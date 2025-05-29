def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = set()
    rectangle_count = 0

    def trace_rectangle(r, c):
        # Find the bottom-right corner of the rectangle
        end_r, end_c = r, c
        while end_c < cols and grid[r][end_c] in ['#', '█']:
            end_c += 1
        while end_r < rows and grid[end_r][c] in ['#', '█']:
            end_r += 1
        # Mark all parts of the rectangle as visited
        for i in range(r, end_r):
            for j in range(c, end_c):
                visited.add((i, j))
        return end_r, end_c

    for r in range(rows):
        for c in range(cols):
            if (r, c) not in visited and grid[r][c] in ['#', '█']:
                trace_rectangle(r, c)
                rectangle_count += 1

    return rectangle_count

# Define the grid
grid = [
    "                                                  ############                  ",
    "                                                  #          #                  ",
    "                                                  #          #                  ",
    "                                                  #          #                  ",
    "                            ######################█##########█###               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "                            #                     #          #  #               ",
    "            ################█#######              #          #  #               ",
    "            #               #######█##############█##########█###               ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              #          #                  ",
    "            #                      #              ############                  ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                            ",
    "            #                      #                                  ##########",
    "            #                      #                                  #        #",
    "            #                      #                                  #        #",
    "            #                      #           #######################█###     #",
    "            #                      #           #                      #  #     #",
    "            #                      #           #                      #  #     #",
    "            #                      #           #######################█###     #",
    "            #                      #                                  #        #",
    "            #                      #                                  ##########",
    "            #                      #                                            ",
    "            #  ####################█####                                        ",
    "            #  #                   #   #                                        ",
    "            ###█####################   #                                        ",
    "               #                       #                                        ",
    "               #                       #                                        ",
    "               #                       #                                        ",
    "               #########################                                        ",
]

# Count the rectangles
rectangle_count = count_rectangles(grid)
print(rectangle_count)