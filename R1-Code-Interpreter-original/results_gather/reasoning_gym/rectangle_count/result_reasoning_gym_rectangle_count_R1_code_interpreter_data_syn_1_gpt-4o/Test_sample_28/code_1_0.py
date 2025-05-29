def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    rectangles = 0

    # Function to check if a rectangle can be formed
    def is_rectangle(r1, c1, r2, c2):
        if r1 >= r2 or c1 >= c2:
            return False
        for r in range(r1, r2 + 1):
            if grid[r][c1] not in '#█' or grid[r][c2] not in '#█':
                return False
        for c in range(c1, c2 + 1):
            if grid[r1][c] not in '#█' or grid[r2][c] not in '#█':
                return False
        return True

    # Check all possible top-left and bottom-right corners
    for r1 in range(rows):
        for c1 in range(cols):
            if grid[r1][c1] in '#█':
                for r2 in range(r1 + 1, rows):
                    for c2 in range(c1 + 1, cols):
                        if grid[r2][c2] in '#█' and is_rectangle(r1, c1, r2, c2):
                            rectangles += 1
                            # If there's an overlap, count the second rectangle
                            if grid[r1][c1] == '█' or grid[r2][c2] == '█':
                                rectangles += 1

    return rectangles

# Define the grid
grid = [
    "                                    #######################################     ",
    "                                    #                                     #     ",
    " ###################################█###################                  #     ",
    " #                                  #                  #                  #     ",
    " #   ##############                 #                  #                  #     ",
    " #   #            #                 #                  #                  #     ",
    " #   #            #                 #                  #                  #     ",
    " #   #            #                 #                  #                  #     ",
    " #   #            #                 #                  #                  #     ",
    " #   #            #                 #                  #                  #     ",
    " #   #            #                 #            ####  #                  #     ",
    " ####█############█#################█############█##█###                  #     ",
    "     #            #                 #            #  #                     #     ",
    "     #            #                 #            #  #                     #     ",
    "     ##############                 #            #  #                     #     ",
    "                                    #            #  #                     #     ",
    "                                    #            #  #                     #     ",
    "                                    #############█##█######################     ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                                                 #  #                           ",
    "                             ####################█##█#########                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             #                   #  #        #                  ",
    "                             ####################█##█#########                  ",
    "                                                 ####                           "
]

# Count rectangles
rectangles_count = count_rectangles(grid)
print(rectangles_count)