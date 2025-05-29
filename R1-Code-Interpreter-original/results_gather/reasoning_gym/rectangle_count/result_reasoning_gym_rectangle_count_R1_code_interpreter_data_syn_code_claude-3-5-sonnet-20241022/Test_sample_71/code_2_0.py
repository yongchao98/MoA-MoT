def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_border(ch):
        return ch in '#█'
    
    def validate_rectangle(x1, y1, x2, y2):
        # Must have minimum size
        if x2 <= x1 + 1 or y2 <= y1 + 1:
            return False
            
        # Check all corners
        if not all(is_border(grid[y][x]) for x, y in [
            (x1, y1), (x2, y1),
            (x1, y2), (x2, y2)
        ]):
            return False
            
        # Check horizontal borders
        for x in range(x1 + 1, x2):
            if not (is_border(grid[y1][x]) and is_border(grid[y2][x])):
                return False
        
        # Check vertical borders
        for y in range(y1 + 1, y2):
            if not (is_border(grid[y][x1]) and is_border(grid[y][x2])):
                return False
        
        # Check if it's a valid rectangle (has some empty space inside)
        has_space = False
        for y in range(y1 + 1, y2):
            for x in range(x1 + 1, x2):
                if grid[y][x] == ' ':
                    has_space = True
                    break
            if has_space:
                break
                
        return has_space
    
    rectangles = set()
    # Find all valid rectangles
    for y1 in range(height-1):
        for x1 in range(width-1):
            if not is_border(grid[y1][x1]):
                continue
            
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if validate_rectangle(x1, y1, x2, y2):
                        rectangles.add((x1, y1, x2, y2))
    
    return len(rectangles)

# Test grid
test_grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                          ###   
                                                                          # #   
                                                                          # #   
                                                                          # #   
                           ################################################ #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                             ## #   
                           #                                  #####      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           #                                  #   #      ## #   
                           ###################################█###█######## #   
                                                              #   #       # #   
                                                              #   #       # #   
                                                              #   #       █#█###
                                                              #   #       █ #  #
                                                              #   #       █ #  #
                                                              #####       █ #  #
                                                                          █ #  #
                                                                          █ #  #
                                                                          █##  #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                                          #    #
                                                #################         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                                                #               #         #    #
                       #########################█████████████████######## #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # #    #
                       #                                                # ######
                       #                       #########                #       
                       #                       #       #                #       
                       #                       #       #                #       
             ##########█████████###############█████████#################       
             #                 #                                                
             ###################                                                """

print(find_rectangles(test_grid))