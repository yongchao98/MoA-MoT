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
]

# Function to count rectangles
def count_rectangles(grid):
    # Initialize rectangle count
    rectangle_count = 0
    in_rectangle = False
    
    # Iterate over each row in the grid
    for row in grid:
        # Check for the start of a rectangle
        if row.strip().startswith('#') and row.strip().endswith('#'):
            if not in_rectangle:
                rectangle_count += 1
                in_rectangle = True
        else:
            in_rectangle = False
        
        # Check for overlapping rectangles
        if '█' in row:
            rectangle_count += 1
    
    return rectangle_count

# Count the rectangles in the grid
total_rectangles = count_rectangles(grid)

# Print the total count of rectangles
print(total_rectangles)