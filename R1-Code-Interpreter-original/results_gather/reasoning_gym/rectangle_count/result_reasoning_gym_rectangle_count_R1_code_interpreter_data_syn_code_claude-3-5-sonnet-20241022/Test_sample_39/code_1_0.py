def find_rectangles(grid):
    if not grid:
        return 0
        
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Function to check if a point is part of a valid rectangle border
    def is_border(y, x):
        return 0 <= y < height and 0 <= x < width and grid[y][x] in ['#', '█']
    
    # Function to validate if coordinates form a valid rectangle
    def is_valid_rectangle(top, left, bottom, right):
        # Check if all four corners exist
        if not all(is_border(y, x) for y, x in [
            (top, left), (top, right), (bottom, left), (bottom, right)
        ]):
            return False
            
        # Check horizontal borders
        for x in range(left, right + 1):
            if not (is_border(top, x) and is_border(bottom, x)):
                return False
                
        # Check vertical borders
        for y in range(top, bottom + 1):
            if not (is_border(y, left) and is_border(y, right)):
                return False
                
        return True

    rectangles = set()
    
    # Find all potential rectangles
    for top in range(height):
        for left in range(width):
            if not is_border(top, left):
                continue
                
            # Skip if not a potential top-left corner
            if (top > 0 and is_border(top-1, left)) or (left > 0 and is_border(top, left-1)):
                continue
                
            # Find matching bottom-right corners
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if is_valid_rectangle(top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

grid = [
"                                                                                ",
"                                                                 ############   ",
"                                                                 #          #   ",
"                                                                 #          #   ",
"                                            #############        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                                            #           #        #          #   ",
"                   #########################█###########█########█########  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        #       #  #   ",
"                   #                        #           #        ########█###   ",
"                   #                        #           #                #      ",
"                   #                        #           #                #      ",
"                   #                        #           #                #      ",
"                   #                        #           #                #      ",
"                   #                        #           #                #      ",
"                   #                        #           #                #      ",
"                   #########################█###########█#################      ",
"                                            #           #                       ",
"                                            #           #                       ",
"                                            #           #                       ",
"                                            #           #                       ",
"                                            #           #                       ",
"                                            #           #                       ",
"                                            #############                       ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                                ",
"                                                                   ########     ",
"                                                                   #      #     ",
"                                                                   #      #     ",
"                                                                   #      #     ",
"                                                                   #      #     ",
"                                                                   #      #     ",
"                                                                   ######## ### ",
"           ############################                                     # # ",
"           #                          #                                     ### ",
"           ############################                                         ",
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
"                                                                                "]

print(find_rectangles(grid))