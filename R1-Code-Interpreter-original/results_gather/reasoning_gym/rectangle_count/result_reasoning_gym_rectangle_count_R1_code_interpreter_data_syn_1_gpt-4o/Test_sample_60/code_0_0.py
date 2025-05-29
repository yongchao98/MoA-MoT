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
    for i, row in enumerate(grid):
        # Find all occurrences of '#' or '█' in the row
        indices = [j for j, char in enumerate(row) if char in '#█']
        
        # If there are at least two indices, it indicates a rectangle
        if len(indices) >= 2:
            # Add the rectangle defined by the start and end indices to the set
            rectangles.add((i, indices[0], indices[-1]))
    
    # Return the number of unique rectangles
    return len(rectangles)

# Count the rectangles in the grid
rectangle_count = count_rectangles(grid)
print(rectangle_count)