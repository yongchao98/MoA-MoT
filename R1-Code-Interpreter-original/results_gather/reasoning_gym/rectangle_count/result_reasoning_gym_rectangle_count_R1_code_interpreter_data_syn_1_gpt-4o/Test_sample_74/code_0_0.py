def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = set()

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '#':
                # Check for top-left corner of a rectangle
                for br in range(r + 1, rows):
                    for bc in range(c + 1, cols):
                        if grid[br][bc] == '#' and grid[r][bc] == '#' and grid[br][c] == '#':
                            # Check if the rectangle is valid
                            is_rectangle = True
                            for i in range(r, br + 1):
                                if grid[i][c] != '#' and grid[i][c] != '█':
                                    is_rectangle = False
                                    break
                            for i in range(c, bc + 1):
                                if grid[r][i] != '#' and grid[r][i] != '█':
                                    is_rectangle = False
                                    break
                            for i in range(r, br + 1):
                                if grid[i][bc] != '#' and grid[i][bc] != '█':
                                    is_rectangle = False
                                    break
                            for i in range(c, bc + 1):
                                if grid[br][i] != '#' and grid[br][i] != '█':
                                    is_rectangle = False
                                    break
                            if is_rectangle:
                                rectangles.add((r, c, br, bc))

    return len(rectangles)

# Define the grid
grid = [
    "                                    ",
    "                                    ",
    "                                    ",
    "                                    ",
    "                                    ",
    "                                    ",
    "                                    ",
    "       ####################################",
    "       #                                  #",
    "       #         #########################█####",
    "       #         #                        #   #",
    "       #         #                        #   #",
    "       #         #                        #   #",
    "       #         #                        #   #",
    "       #         #                        #   #",
    "       #         #########################█####",
    "       #                                  #",
    "       #                                  #",
    "       #                                  #",
    "       #                                  #",
    "       #                                  #",
    "       #                                  #",
    "       #                                  #",
    "       #                                  #                  #######",
    "       #                                  #                  #     #",
    "       #                                  #                  #     #",
    "       #                                  #                  #     #",
    "       #                                  #                  #     #",
    "       #                                  #                  #     #",
    "       #                                  #                  #     #",
    "       #                                  #                  #     #",
    "       #                                  #                  #     #",
    "       #                                  #                  #     #",
    "       ####################################                  #     #",
    "                                                             #     #",
    "                                                             #     #",
    "                                                             #     #",
    "                                                             #     #",
    "                                                             #     #",
    "                                                             #     #",
    "                                                             #     #",
    "                                    ################################",
    "                                    #                       #",
    "                                    #         ##############█###################",
    "                                    #         #             #                  #",
    "                                    #         #             #                  #",
    "                                    #         #             #                  #",
    "                                    #         #             #                  #",
    "                                    #         #             #                  #",
    "                                    #         ##############█###################",
    "                                    #                       #",
    "                                    #                       #",
    "                                    #                       #",
    "                                    #                       #",
    "                                    #########################  ########",
    "                                                               #      #",
    "                                                               #      #",
    "                                                               #      #",
    "                                    #################          #      #",
    "                                    #               #          #      #",
    "                                    #               #          #      #",
    "                                    #               #          #      #",
    "                                    #               #          #      #",
    "                                    #               #          #      #",
    "                                    #               #          #      #",
    "                                    #               #          #      #",
    "  ##################################█###############█##########█######█#####",
    "  #                                 #               #          #      #    #",
    "  #                                 #################          #      #    #",
    "  #                                                            #      #    #",
    "  #                                                            ########    #",
    "  #                                                                        #",
    "  #                                                                        #",
    "  #                                                                        #",
    "  #                                                                        #",
    "  ##########################################################################",
    "                                    ",
    "                                    ",
    "                                    ",
]

# Count the rectangles
print(count_rectangles(grid))