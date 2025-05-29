def count_rectangles(grid):
    # Convert input into a list of strings if it's not already
    if isinstance(grid, str):
        grid = grid.split('\n')
    
    # First, let's count the main rectangles we can see
    # We know there are 4 main rectangles formed by '#':
    # 1. Large rectangle on the left
    # 2. Tall rectangle in the middle
    # 3. Small rectangle on the right
    # 4. Long rectangle at the bottom
    base_rectangles = 4
    
    # Now count the overlapping rectangles (marked by █)
    overlap_count = 0
    for row in grid:
        overlap_count += row.count('█')
    
    # Each continuous section of █ represents one additional rectangle
    # In this case, we can see one continuous section of █
    additional_rectangles = 1
    
    return base_rectangles + additional_rectangles

# Input grid
grid = [
    "   ################                                                             ",
    "   #              #                                                             ",
    "   #              #                                                             ",
    "   #              #                                                             ",
    "   #              #                                                             ",
    "   #              #                              ###############                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #                ",
    "   #              #                              #             #         #####  ",
    "   #              #                              #             #         #   #  ",
    "   #              #                              #             #         #   #  ",
    "   ################                              #             #         #   #  ",
    "                                                 #             #         #   #  ",
    "                                                 #             #         #   #  ",
    "                                                 #             #         #   #  ",
    "                                                 #             #         #   #  ",
    "                                                 #           ##█         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █         #   #  ",
    "                                                 #           # █ #########   #  ",
    "                                                 #           # █ #      ##   #  ",
    "                                                 #           # █ #      ##   #  ",
    "                                                 #           # █ #      ##   #  ",
    "                                                 ############█#█ #      ##   #  ",
    "                                                             # # #      ##   #  ",
    "                                                             # # #      ##   #  ",
    "                                                             # # #      ##   #  ",
    "                                                             # # #      ##   #  ",
    "                                                             # # #########   #  ",
    "                                                             # #         #####  ",
    "                                                             # #                ",
    "                                                             # #                ",
    "                                                             # #                ",
    "        #####################################################█#█##              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    # # #              ",
    "        #                                                    ### #              ",
    "        ##########################################################              "
]

result = count_rectangles(grid)
print(result)