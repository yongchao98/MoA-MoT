def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    visited = set()
    rectangles = set()
    
    def is_border(y, x):
        if 0 <= y < height and 0 <= x < width:
            return grid[y][x] in '#█'
        return False
    
    def find_rectangle(start_y, start_x):
        if not is_border(start_y, start_x):
            return None
            
        # Find right edge
        end_x = start_x
        while end_x < width and is_border(start_y, end_x):
            end_x += 1
        end_x -= 1
        
        # Find bottom edge
        end_y = start_y
        while end_y < height and is_border(end_y, start_x):
            end_y += 1
        end_y -= 1
        
        # Validate rectangle
        for y in range(start_y, end_y + 1):
            if not (is_border(y, start_x) and is_border(y, end_x)):
                return None
                
        for x in range(start_x, end_x + 1):
            if not (is_border(start_y, x) and is_border(end_y, x)):
                return None
                
        return (start_y, start_x, end_y, end_x)
    
    # Find all rectangles
    for y in range(height):
        for x in range(width):
            if (y, x) not in visited and is_border(y, x):
                # Check if it's a top-left corner
                if (not is_border(y-1, x)) and (not is_border(y, x-1)):
                    rect = find_rectangle(y, x)
                    if rect:
                        start_y, start_x, end_y, end_x = rect
                        # Mark all border points as visited
                        for ry in range(start_y, end_y + 1):
                            for rx in range(start_x, end_x + 1):
                                if is_border(ry, rx):
                                    visited.add((ry, rx))
                        rectangles.add(rect)
    
    return len(rectangles)

# Input grid
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