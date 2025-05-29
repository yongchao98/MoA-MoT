def count_rectangles(grid):
    rectangles = 0
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '#':
                # Check if this is the top-left corner of a rectangle
                if (r == 0 or grid[r-1][c] != '#') and (c == 0 or grid[r][c-1] != '#'):
                    rectangles += 1
    return rectangles

# Define the grid
grid = [
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "        ###############################################                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        #                                             #                         ",
    "        ###############################################                         ",
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
    "                              ##################                                ",
    "                              #                #                                ",
    "                              #                #                                ",
    "                              #                #                                ",
    "                              ##################                                ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                                                 ",
    "                                                        ##############           ",
    "                                                        #            #           ",
    "                                                        ##############           ",
    "                                                                                 "
]

# Count the rectangles
rectangles_count = count_rectangles(grid)
print(rectangles_count)