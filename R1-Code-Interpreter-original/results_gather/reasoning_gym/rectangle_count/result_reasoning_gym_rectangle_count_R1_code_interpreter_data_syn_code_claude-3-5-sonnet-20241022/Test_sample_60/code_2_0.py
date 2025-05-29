def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_valid_char(c):
        return c in '#█'
    
    def find_rectangle_end(start_y, start_x):
        # Find right end of top line
        right = start_x
        while right < width and is_valid_char(grid[start_y][right]):
            right += 1
        right -= 1
        
        if right <= start_x:
            return None
        
        # Find bottom end
        bottom = start_y
        while bottom < height and is_valid_char(grid[bottom][start_x]):
            bottom += 1
        bottom -= 1
        
        if bottom <= start_y:
            return None
            
        # Verify rectangle
        # Check right vertical line
        for y in range(start_y, bottom + 1):
            if not is_valid_char(grid[y][right]):
                return None
                
        # Check bottom horizontal line
        for x in range(start_x, right + 1):
            if not is_valid_char(grid[bottom][x]):
                return None
                
        # Check if interior is valid (space or overlap markers)
        for y in range(start_y + 1, bottom):
            for x in range(start_x + 1, right):
                if grid[y][x] not in ' #█':
                    return None
        
        return (bottom, right)
    
    rectangles = set()
    # Only look for top-left corners that start with '#' or '█'
    for y in range(height):
        for x in range(width):
            if is_valid_char(grid[y][x]):
                # Check if this could be a top-left corner
                if (y == 0 or not is_valid_char(grid[y-1][x])) and \
                   (x == 0 or not is_valid_char(grid[y][x-1])):
                    result = find_rectangle_end(y, x)
                    if result:
                        bottom, right = result
                        rectangles.add((y, x, bottom, right))
    
    return len(rectangles)

# Input grid (using the same grid as before)
grid_str = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
         ##############################################                         
         #                                            #                         
         #                                            #                         
         #                                            #                         
         #                                            #                         
         #                                            #                         
         #                                            #                         
         ##############################################                         
                     ################################################           
                     #                                              #           
                     #                                              #           
                     #                                              #           
                     #                                              #           
            #########█##############################################█######     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #  ###
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #      ##█##############################################█###  #  # #
            #      # ################################################  #  #  # #
            #      #                                                   #  #  # #
            #      #                                                   #  #  # #
            #      #####################################################  #  # #
            #                                                             #  # #
            #                                                             #  # #
            #                                                             #  # #
            #                                                             #  # #
            #                                                             #  # #
            #                                                             #  ###
            #                                                     ########█###  
            #                                                     #       #  #  
            #                                                     #       #  #  
            #                                                     #       #  #  
   #########█##############                                       #       #  #  
   #        ##############█#######################################█########  #  
   #                      #                                       ############  
   ########################                                                     
                                                                                
                                                                                
"""

# Clean up the grid
grid = [list(line) for line in grid_str.splitlines() if line.strip()]

print(find_rectangles(grid))