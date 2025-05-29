def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = set()

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in '#█':
                # Try to find a rectangle starting at (r, c)
                for rr in range(r + 1, rows):
                    if grid[rr][c] not in '#█':
                        break
                    for cc in range(c + 1, cols):
                        if grid[r][cc] not in '#█':
                            break
                        if grid[rr][cc] in '#█':
                            # Check if we have a complete rectangle
                            if all(grid[rr][x] in '#█' for x in range(c, cc + 1)) and all(grid[x][cc] in '#█' for x in range(r, rr + 1)):
                                rectangles.add((r, c, rr, cc))
    return len(rectangles)

# Define the grid
grid = [
    "                                                                     ########## ",
    "                                                                     #        # ",
    "                     ######################################          #        # ",
    "                     #                          ##########█########  #        # ",
    "                     #                          #         #       #  #        # ",
    "                     #                          #         #       #  #        # ",
    "                     ###########################█##########       #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                        ########################█#################█##█######  # ",
    "                        #                       #                 #  #     #  # ",
    "                        ########################█#################█##█######  # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #                 #  #        # ",
    "                                                #      ###########█##█###     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  #  #     # ",
    "                                                #      #          #  ###█###### ",
    "                                                #      #          #     #       ",
    "                                                #      #          #     #       ",
    "                                                #      #          #     #       ",
    "                                                #######█###########     #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                                       #                #       ",
    "                                          #############█################█       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "                                          #            #                █       ",
    "       ##############                     #############█########████####█       ",
    "       #            #                                  #########█##█#####       ",
    "       ##############                                           ####            "
]

# Count the rectangles
print(count_rectangles(grid))