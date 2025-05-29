def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = set()

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in '#█':
                # Check if this is a top-left corner of a rectangle
                if (r == 0 or grid[r-1][c] not in '#█') and (c == 0 or grid[r][c-1] not in '#█'):
                    # Find the bottom-right corner of the rectangle
                    rr, cc = r, c
                    while rr + 1 < rows and grid[rr + 1][c] in '#█':
                        rr += 1
                    while cc + 1 < cols and grid[r][cc + 1] in '#█':
                        cc += 1
                    # Verify the rectangle
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