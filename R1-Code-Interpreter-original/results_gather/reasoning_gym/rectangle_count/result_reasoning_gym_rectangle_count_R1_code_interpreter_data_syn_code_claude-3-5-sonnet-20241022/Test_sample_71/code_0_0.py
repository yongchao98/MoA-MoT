def find_rectangles(grid):
    if not grid:
        return 0
    
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def is_corner(x, y):
        return grid[y][x] in '#█'
    
    def verify_rectangle(x1, y1, x2, y2):
        # Check if this forms a valid rectangle
        if x2 <= x1 or y2 <= y1:
            return False
            
        # Check horizontal lines
        for x in range(x1, x2 + 1):
            if not (grid[y1][x] in '#█' and grid[y2][x] in '#█'):
                return False
                
        # Check vertical lines
        for y in range(y1, y2 + 1):
            if not (grid[y][x1] in '#█' and grid[y][x2] in '#█'):
                return False
                
        return True
    
    rectangles = set()
    # Find all possible rectangles
    for y1 in range(height):
        for x1 in range(width):
            if not is_corner(x1, y1):
                continue
                
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if is_corner(x2, y2) and verify_rectangle(x1, y1, x2, y2):
                        rectangles.add((x1, y1, x2, y2))
    
    return len(rectangles)

# Test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
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

print(find_rectangles(grid))