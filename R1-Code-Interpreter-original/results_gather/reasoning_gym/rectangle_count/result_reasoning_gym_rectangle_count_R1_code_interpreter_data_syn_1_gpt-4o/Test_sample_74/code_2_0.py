def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = set()

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '#':
                # Check for top-left corner of a rectangle
                for br in range(r + 1, rows):
                    if grid[br][c] not in ['#', '█']:
                        break
                    for bc in range(c + 1, cols):
                        if grid[br][bc] not in ['#', '█']:
                            break
                        if grid[r][bc] not in ['#', '█']:
                            continue
                        # Check if the rectangle is valid
                        is_rectangle = True
                        for i in range(r, br + 1):
                            if grid[i][c] not in ['#', '█']:
                                is_rectangle = False
                                break
                        for i in range(c, bc + 1):
                            if grid[r][i] not in ['#', '█']:
                                is_rectangle = False
                                break
                        for i in range(r, br + 1):
                            if grid[i][bc] not in ['#', '█']:
                                is_rectangle = False
                                break
                        for i in range(c, bc + 1):
                            if grid[br][i] not in ['#', '█']:
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