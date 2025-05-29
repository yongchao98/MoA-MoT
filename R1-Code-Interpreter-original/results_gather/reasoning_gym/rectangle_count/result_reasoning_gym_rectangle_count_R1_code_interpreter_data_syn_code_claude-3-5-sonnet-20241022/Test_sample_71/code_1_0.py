def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line]
    
    # Create a 2D array of the grid
    grid_array = []
    for line in lines:
        grid_array.append(list(line))
    
    height = len(grid_array)
    width = len(grid_array[0])
    
    def is_valid_rect(x1, y1, x2, y2):
        # Check minimum size
        if x2 <= x1 + 1 or y2 <= y1 + 1:
            return False
            
        # Check if corners are borders
        corners = [(x1, y1), (x2, y1), (x1, y2), (x2, y2)]
        for x, y in corners:
            if grid_array[y][x] not in '#█':
                return False
        
        # Check horizontal borders
        for x in range(x1, x2 + 1):
            if grid_array[y1][x] not in '#█' or grid_array[y2][x] not in '#█':
                return False
        
        # Check vertical borders
        for y in range(y1, y2 + 1):
            if grid_array[y][x1] not in '#█' or grid_array[y][x2] not in '#█':
                return False
        
        # Check if rectangle has interior space
        has_space = False
        for y in range(y1 + 1, y2):
            for x in range(x1 + 1, x2):
                if grid_array[y][x] == ' ':
                    has_space = True
                    break
            if has_space:
                break
        
        return has_space
    
    rectangles = set()
    # Find all possible rectangles
    for y1 in range(height - 1):
        for x1 in range(width - 1):
            if grid_array[y1][x1] not in '#█':
                continue
            
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if grid_array[y2][x2] not in '#█':
                        continue
                    
                    if is_valid_rect(x1, y1, x2, y2):
                        # Check if this rectangle is not contained within another
                        is_unique = True
                        for rx1, ry1, rx2, ry2 in rectangles:
                            if (rx1 <= x1 and x2 <= rx2 and 
                                ry1 <= y1 and y2 <= ry2 and 
                                (rx1 != x1 or rx2 != x2 or ry1 != y1 or ry2 != y2)):
                                is_unique = False
                                break
                        if is_unique:
                            rectangles.add((x1, y1, x2, y2))
    
    return len(rectangles)

# Your grid
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