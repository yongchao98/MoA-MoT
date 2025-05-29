# Define the grid as a list of strings
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "         ##############################################                         ",
    "         #                                            #                         ",
    "         #                                            #                         ",
    "         #                                            #                         ",
    "         #                                            #                         ",
    "         #                                            #                         ",
    "         #                                            #                         ",
    "         ##############################################                         ",
    "                     ################################################           ",
    "                     #                                              #           ",
    "                     #                                              #           ",
    "                     #                                              #           ",
    "                     #                                              #           ",
    "            #########█##############################################█######     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #     ",
    "            #        #                                              #     #  ###",
    "            #        #                                              #     #  # #",
    "            #        #                                              #     #  # #",
    "            #        #                                              #     #  # #",
    "            #        #                                              #     #  # #",
    "            #        #                                              #     #  # #",
    "            #        #                                              #     #  # #",
    "            #      ##█##############################################█###  #  # #",
    "            #      # ################################################  #  #  # #",
    "            #      #                                                   #  #  # #",
    "            #      #                                                   #  #  # #",
    "            #      #####################################################  #  # #",
    "            #                                                             #  # #",
    "            #                                                             #  # #",
    "            #                                                             #  # #",
    "            #                                                             #  # #",
    "            #                                                             #  # #",
    "            #                                                             #  ###",
    "            #                                                     ########█###  ",
    "            #                                                     #       #  #  ",
    "            #                                                     #       #  #  ",
    "            #                                                     #       #  #  ",
    "   #########█##############                                       #       #  #  ",
    "   #        ##############█#######################################█########  #  ",
    "   #                      #                                       ############  ",
    "   ########################                                                      ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
]

# Function to count rectangles
def count_rectangles(grid):
    # Initialize a set to store unique rectangles
    rectangles = set()
    
    # Iterate over each row in the grid
    for i in range(len(grid) - 1):
        for j in range(len(grid[i]) - 1):
            # Check for the top-left corner of a rectangle
            if grid[i][j] in '#█' and grid[i][j+1] in '#█' and grid[i+1][j] in '#█':
                # Find the bottom-right corner of the rectangle
                for k in range(i, len(grid)):
                    if grid[k][j] not in '#█':
                        break
                    for l in range(j, len(grid[k])):
                        if grid[k][l] not in '#█':
                            break
                        # Check if this is a bottom-right corner
                        if grid[k][j] in '#█' and grid[i][l] in '#█':
                            rectangles.add((i, j, k, l))
    
    # Return the number of unique rectangles
    return len(rectangles)

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)