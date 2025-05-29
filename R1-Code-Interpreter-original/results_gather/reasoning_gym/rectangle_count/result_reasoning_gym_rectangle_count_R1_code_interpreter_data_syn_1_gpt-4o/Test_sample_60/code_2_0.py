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
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            # Check for the top-left corner of a rectangle
            if grid[i][j] == '#':
                # Determine the width of the rectangle
                width = 0
                while j + width < len(grid[i]) and grid[i][j + width] in '#█':
                    width += 1
                
                # Determine the height of the rectangle
                height = 0
                while i + height < len(grid) and grid[i + height][j] in '#█':
                    height += 1
                
                # Add the rectangle defined by the top-left corner and its dimensions to the set
                rectangles.add((i, j, width, height))
    
    # Return the number of unique rectangles
    return len(rectangles)

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)