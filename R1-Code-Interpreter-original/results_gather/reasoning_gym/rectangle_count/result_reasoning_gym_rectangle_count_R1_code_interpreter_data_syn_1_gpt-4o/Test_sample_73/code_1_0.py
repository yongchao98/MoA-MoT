def count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0]) if rows > 0 else 0
    visited = set()
    rectangles = 0

    # Function to check if a rectangle is formed
    def is_rectangle(r1, c1, r2, c2):
        if r1 >= r2 or c1 >= c2:
            return False
        # Check top and bottom edges
        for c in range(c1, c2 + 1):
            if grid[r1][c] not in '#█' or grid[r2][c] not in '#█':
                return False
        # Check left and right edges
        for r in range(r1, r2 + 1):
            if grid[r][c1] not in '#█' or grid[r][c2] not in '#█':
                return False
        return True

    # Iterate over all possible top-left corners
    for r1 in range(rows):
        for c1 in range(cols):
            if grid[r1][c1] in '#█' and (r1, c1) not in visited:
                # Iterate over all possible bottom-right corners
                for r2 in range(r1 + 1, rows):
                    for c2 in range(c1 + 1, cols):
                        if grid[r2][c2] in '#█' and is_rectangle(r1, c1, r2, c2):
                            # Mark all positions of this rectangle as visited
                            for r in range(r1, r2 + 1):
                                for c in range(c1, c2 + 1):
                                    visited.add((r, c))
                            rectangles += 1

    return rectangles

# Define the grid
grid = [
    "                                            #########                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #                           ",
    "                                            #       #             ##############",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                       #####█#######█#############█########    #",
    "                                       #    #       #             #       #    #",
    "                                       #    #       #             #       #    #",
    "                                       #    #       #             #       #    #",
    "                                       #    #       #             #       #    #",
    "                                       #####█#######█#############█########    #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #       #             #            #",
    "                                            #########             #            #",
    "                                                                  ##############",
    "                                                                  ##          ##",
    "                    ##############################################██######### ##",
    "                    #                    ##########               ##        # ##",
    "                    #                    #        #               ##        # ##",
    "                    #                    #        #               ##        # ##",
    "                    #                    #        #               ##        # ##",
    "                    #                    #        #               ##        # ##",
    "                  ##█######              #        #               ##        # ##",
    "                  # #     #              #        #               ##        # ##",
    "                  # #     #              #        #               ##        # ##",
    "                  # #     #              #        #               #█########█#█#",
    "                  # #     #              #        #                #        # # ",
    "                  # ######█##############█########█################█########█## ",
    "                  # #█####█##############█########█################█########### ",
    "                  #  #    #              #        #                #         ## ",
    "                  #  #    #              ##########                #         ## ",
    "                  #  #    #                                        #         ## ",
    "                  #  #    #                                        #         ## ",
    "                  #  #    #                                        #         ## ",
    "                  #  #    #                                        ##########█# ",
    "                  ###█#####                                                  #  ",
    "                     #                            ####################       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            #                  #       #  ",
    "                     #                            ############████████########  ",
    "                     #########################################█#############█#  ",
    "                                                              #             #   ",
    "                                                              ###############   ",
    "                                                                                 "
]

# Count the rectangles
rectangles_count = count_rectangles(grid)
print(rectangles_count)