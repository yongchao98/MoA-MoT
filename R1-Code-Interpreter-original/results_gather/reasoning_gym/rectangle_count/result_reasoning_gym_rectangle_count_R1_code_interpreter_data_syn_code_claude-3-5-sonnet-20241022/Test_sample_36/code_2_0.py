def find_rectangles(grid):
    if not grid:
        return 0
        
    height = len(grid)
    width = len(grid[0])
    
    def is_border(y, x):
        if 0 <= y < height and 0 <= x < width:
            return grid[y][x] in '#█'
        return False
    
    def check_rectangle(y1, x1):
        # Find the right edge
        x2 = x1
        while x2 < width and is_border(y1, x2):
            x2 += 1
        x2 -= 1
        
        if x2 <= x1:
            return None
            
        # Find the bottom edge
        y2 = y1
        while y2 < height and is_border(y2, x1):
            y2 += 1
        y2 -= 1
        
        if y2 <= y1:
            return None
            
        # Verify all borders exist
        for x in range(x1, x2 + 1):
            if not (is_border(y1, x) and is_border(y2, x)):
                return None
                
        for y in range(y1, y2 + 1):
            if not (is_border(y, x1) and is_border(y, x2)):
                return None
                
        # Check if it's a valid rectangle (no other borders inside)
        for y in range(y1 + 1, y2):
            for x in range(x1 + 1, x2):
                if is_border(y, x) and not (y in [y1, y2] or x in [x1, x2]):
                    return None
                    
        return (y1, x1, y2, x2)
    
    rectangles = set()
    # Look for top-left corners
    for y in range(height):
        for x in range(width):
            if is_border(y, x):
                # Check if it's a potential top-left corner
                if (not is_border(y-1, x) or y == 0) and (not is_border(y, x-1) or x == 0):
                    rect = check_rectangle(y, x)
                    if rect:
                        rectangles.add(rect)
    
    return len(rectangles)

# Input grid as a list of strings
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                     ##########",
    "                                                                     #        #",
    "                     ######################################          #        #",
    "                     #                          ##########█########  #        #",
    "                     #                          #         #       #  #        #",
    "                     #                          #         #       #  #        #",
    "                     ###########################█##########       #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                        ########################█#################█##█######  #",
    "                        #                       #                 #  #     #  #",
    "                        ########################█#################█##█######  #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #                 #  #        #",
    "                                                #      ###########█##█###     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  #  #     #",
    "                                                #      #          #  ###█######",
    "                                                #      #          #     #      ",
    "                                                #      #          #     #      ",
    "                                                #      #          #     #      ",
    "                                                #######█###########     #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                                       #                #      ",
    "                                          #############█################█      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "                                          #            #                █      ",
    "       ##############                     #############█########████####█      ",
    "       #            #                                  #########█##█#####      ",
    "       ##############                                           ####           "
]

print(find_rectangles(grid))