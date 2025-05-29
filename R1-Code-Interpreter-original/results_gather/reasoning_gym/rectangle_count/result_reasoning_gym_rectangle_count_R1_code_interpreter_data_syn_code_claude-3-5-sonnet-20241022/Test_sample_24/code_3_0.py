def find_complete_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    def is_valid_line_horizontal(y, x1, x2):
        for x in range(x1, x2 + 1):
            if grid[y][x] not in '#█':
                return False
        return True
    
    def is_valid_line_vertical(x, y1, y2):
        for y in range(y1, y2 + 1):
            if grid[y][x] not in '#█':
                return False
        return True
    
    def is_valid_rectangle(top, left, bottom, right):
        # Check all four sides
        if not is_valid_line_horizontal(top, left, right):  # Top
            return False
        if not is_valid_line_horizontal(bottom, left, right):  # Bottom
            return False
        if not is_valid_line_vertical(left, top, bottom):  # Left
            return False
        if not is_valid_line_vertical(right, top, bottom):  # Right
            return False
        return True
    
    # Find potential corners
    corners = []
    for i in range(height):
        for j in range(width):
            if grid[i][j] in '#█':
                corners.append((i, j))
    
    # Check each possible rectangle
    for i, (top, left) in enumerate(corners):
        for bottom, right in corners[i:]:
            if bottom > top and right > left:  # Ensure correct orientation
                if is_valid_rectangle(top, left, bottom, right):
                    # Verify it's a complete rectangle
                    rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

# Input grid as a list of strings
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                 ###########################                                    ",
    "                 #                         #                                    ",
    "                 #                         #                                    ",
    "                 #                         #                             ####   ",
    "                 #                         #                             #  #   ",
    "                 #                         #   ######################    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 #                         #   #                    #    #  #   ",
    "                 ###########################   #                    #    #  #   ",
    "                                               #                    #    #  #   ",
    "                                ###############█####################█##  #  #   ",
    "                                #              #                    # #  #  #   ",
    "                                #              #                    # #  #  #   ",
    "                                ###############█####################███##█##█#  ",
    "                                               #                    ##   #  ##  ",
    "                                               # ###################██###█##██##",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  ##   #  ## #",
    "                                               # #                  #####█##█# #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # #                  #    #  #  #",
    "                                               # ###################█####█##█###",
    "                                               #                    #    #  #   ",
    "                                               #                    #    #  #   ",
    "                                               #                    #    #  #   ",
    "                                               #                    #    #  #   ",
    "                                               #        ########    #    #  #   ",
    "                                               #        #      #    #    #  #   ",
    "                                               #        #      #    #    #  #   ",
    "                                               #        #      #    #    #  #   ",
    "                                               #########█######█#####    #  #   ",
    "                                                        #      #         #  #   ",
    "                                                        ########         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    " ###############################################################         #  #   ",
    " #                                                             #         #  #   ",
    " #                                                             #         #  #   ",
    " #                                                             #         #  #   ",
    " #                                                             #         #  #   ",
    " #          ############################                       #         #  #   ",
    " ###########█##########################█########################         #  #   ",
    "            #                          #                                 #  #   ",
    "            ############################                                 #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         #  #   ",
    "                                                                         ####   "
]

print(find_complete_rectangles(grid))