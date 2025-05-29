def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = 0

    # Function to check if a position is a top-left corner of a rectangle
    def is_top_left(r, c):
        if grid[r][c] not in ('#', '█'):
            return False
        if r > 0 and grid[r-1][c] in ('#', '█'):
            return False
        if c > 0 and grid[r][c-1] in ('#', '█'):
            return False
        return True

    # Function to check if a position is a bottom-right corner of a rectangle
    def is_bottom_right(r, c):
        if grid[r][c] not in ('#', '█'):
            return False
        if r < rows - 1 and grid[r+1][c] in ('#', '█'):
            return False
        if c < cols - 1 and grid[r][c+1] in ('#', '█'):
            return False
        return True

    # Scan the grid to count rectangles
    for r in range(rows):
        for c in range(cols):
            if is_top_left(r, c):
                # Find the bottom-right corner of the rectangle
                for br in range(r, rows):
                    for bc in range(c, cols):
                        if is_bottom_right(br, bc):
                            rectangles += 1
                            break
                    else:
                        continue
                    break

    return rectangles

# Define the grid
grid = [
    "                                                 ##############                 ",
    "                                                 #            #                 ",
    "                                                 #            #                 ",
    "                                                 #            #                 ",
    "                                       ##########█############█                 ",
    "                                       #         #            █                 ",
    "                                       #         #            █                 ",
    "                                       #         #            █                 ",
    "                                       #         #            █                 ",
    "                                       #         #            █                 ",
    "                                       #         #            █                 ",
    "                                       #         #            █                 ",
    "       ################################█#########█#########   █                 ",
    "       #                               #         #        #   █                 ",
    "       #                               #         #        #   █                 ",
    "       #                               #         #        #   █                 ",
    "       #                               #         #        #   █           ##### ",
    "       #                               #         #        #   █           #   # ",
    "       #                               #         #        #   █           #   # ",
    "       #                               #         #        #   █           #   # ",
    "       #                               #     #####        #   █           #   # ",
    "       #                               #     #  ##        #   █           #   # ",
    "       #                               #     #  ##        #   █           #   # ",
    "       █###############################█#####█████####### #   █           #   # ",
    "       █                               #         #      # #   █           #   # ",
    "       █                               #         #      # #   █           #   # ",
    "       █           #########           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       #           #         #      # #   █           #   # ",
    "       █           #       █###########█###      #      # #   █           #   # ",
    "       █           #       █           #  #      #      # #   █           #   # ",
    "       █           ########█           #  #      #      # #   █           #   # ",
    "       █###################█###########█##█######█######█##   █           #   # ",
    "       #                   #           #  #      #      #     █           #   # ",
    "       #                   #           #  #      #      #     █           #   # ",
    "       #                   #           #  #      #      #     █           #   # ",
    "       #                   #           #  #      #######█#####█           #   # ",
    "       #                   #           #  #             #     #           #   # ",
    "       #                   #           #  #             #     #           #   # ",
    "#######█##############     #           #  #             #     #           #   # ",
    "#      #             #     #           #  #             #     #           #   # ",
    "#      ##############█#####█###########█##█##############     #           #   # ",
    "######################     #           #  #                   #           #   # ",
    "                           #           #  #                   #           #   # ",
    "                           #           #  #                   #           #   # ",
    "                           #           #  #                   #           #   # ",
    "                           #           #  #                   #           #   # ",
    "                           #           #  #                   #           #   # ",
    "                           ############█###                   #           #   # ",
    "                                       #                      #           ##### ",
    "                                       #                      #                 ",
    "                                       ########################                 "
]

# Count the rectangles
rectangles_count = count_rectangles(grid)
print(rectangles_count)