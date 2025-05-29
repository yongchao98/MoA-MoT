def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0])
    visited = set()
    rectangles = set()

    def is_valid_rectangle(r1, c1, r2, c2):
        # Check top and bottom borders
        for c in range(c1, c2 + 1):
            if grid[r1][c] not in ('#', '█') or grid[r2][c] not in ('#', '█'):
                return False
        # Check left and right borders
        for r in range(r1, r2 + 1):
            if grid[r][c1] not in ('#', '█') or grid[r][c2] not in ('#', '█'):
                return False
        return True

    for r in range(rows):
        for c in range(cols):
            if grid[r][c] in ('#', '█') and (r, c) not in visited:
                # Find the bottom-right corner
                br_r, br_c = r, c
                while br_r + 1 < rows and grid[br_r + 1][c] in ('#', '█'):
                    br_r += 1
                while br_c + 1 < cols and grid[r][br_c + 1] in ('#', '█'):
                    br_c += 1
                # Check if it's a valid rectangle
                if is_valid_rectangle(r, c, br_r, br_c):
                    rectangles.add((r, c, br_r, br_c))
                    # Mark all cells in this rectangle as visited
                    for i in range(r, br_r + 1):
                        for j in range(c, br_c + 1):
                            visited.add((i, j))

    return len(rectangles)

# Define the grid
grid = [
    "                                   ###################################          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "                                   #                                 #          ",
    "  #####                            #                                 #          ",
    "  #   #                            #                                 #          ",
    "  #   #                            #                                 #          ",
    "  #   #                            #                                 #          ",
    "  #   #                            #                                 #          ",
    "  #   #                            #                                 #          ",
    "  #   #                            #                                 #          ",
    "  #   #                            #                                 #######    ",
    "  #   #                            █##############################   ##    #    ",
    "  #   #                            █                             #   ##    #    ",
    "  #   #                            █                             #   ##    #    ",
    "  #   #                            █                             #   ##    #    ",
    "  #   #                            █                             #   ##    #    ",
    "  #   #                            █                             #   ##    #    ",
    "  #   #                            █   ##############            #   ##    #    ",
    "  #   #                            █   #            #            #   ##    #    ",
    "  #  #█#################           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            #   ##    #    ",
    "  #  ##                #           █   #            #            █#####    #    ",
    "  #  ##                #           █   #            #            █  ###    #    ",
    "  #  ##                #           █   #            #            █  ###    #    ",
    "  #  ##                #           █###█############█############█  ###    #    ",
    "  #  ##                #           #   #            #            #  ###    #    ",
    "  #  ##                #           #   #            #            #  ###    #    ",
    "  #  ##                #           #   #            #            #  ###    #    ",
    "  #  ##                #           #   #            #            #  ###    #    ",
    "  #  ##                #           #   #            #            #  ###    #    ",
    "  #  ##                #           #   #            #            ######    #    ",
    "  #  ##                #           #   #            #                ##    #    ",
    "  #  ##                #           #   ##############                ##    #    ",
    "  #  ##                #           #                                 ##    #    ",
    "  #  ##            ####█###########█########################         ##    #    ",
    "  #  ##            #   #           #                       #         ##    #    ",
    "  #  ##            #   #          #█#################      #         ##    #    ",
    "  #  ##            #   #          ##                #      #         ##    #    ",
    "  #  ##            #   #          #█#################      #         ##    #    ",
    "  #  ##        ####█###█###########█##########             #         #######    ",
    "  #  #█########█###█####           #         #             #         #          ",
    "  #   #        #   #               #         #             #         #          ",
    "  #   #        #   #               #         #             #         #          ",
    "  #   #        #   ################█#########█##############         #          ",
    "  #   #        #                   #         #                       #          ",
    "  #   #        #                   #         #                       #          ",
    "  #####        #                   #         #                       #          ",
    "               #                   #         #                       #          ",
    "               #                   ##########█########################          ",
    "               #                             #                                  ",
    "               #                             #                                  ",
    "               #                             #                                  ",
    "               #                             #                                  ",
    "               #                             #                                  ",
    "               #                             #                                  ",
    "               #                             #                                  ",
    "               #                             #                                  ",
    "               ###############################                                  ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                "
]

# Count the rectangles
print(count_rectangles(grid))