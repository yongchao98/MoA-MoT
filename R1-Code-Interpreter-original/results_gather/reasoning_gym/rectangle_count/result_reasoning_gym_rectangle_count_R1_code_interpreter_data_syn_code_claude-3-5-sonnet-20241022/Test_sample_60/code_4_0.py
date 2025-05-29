def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_border(c):
        return c in '#█'
    
    def find_complete_rectangle(y, x):
        # Skip if not a potential top-left corner
        if not is_border(grid[y][x]):
            return None
        if y > 0 and is_border(grid[y-1][x]):  # Check no border above
            return None
        if x > 0 and is_border(grid[y][x-1]):  # Check no border to left
            return None
            
        # Find right end of top border
        right = x + 1
        while right < width and is_border(grid[y][right]):
            right += 1
        right -= 1
        
        if right == x:  # Not wide enough to be a rectangle
            return None
            
        # Find bottom end of left border
        bottom = y + 1
        while bottom < height and is_border(grid[bottom][x]):
            bottom += 1
        bottom -= 1
        
        if bottom == y:  # Not tall enough to be a rectangle
            return None
            
        # Verify right border
        for y_pos in range(y, bottom + 1):
            if not is_border(grid[y_pos][right]):
                return None
                
        # Verify bottom border
        for x_pos in range(x, right + 1):
            if not is_border(grid[bottom][x_pos]):
                return None
                
        # Verify interior (must be space or overlap markers)
        for y_pos in range(y + 1, bottom):
            for x_pos in range(x + 1, right):
                if grid[y_pos][x_pos] not in ' #█':
                    return None
                    
        return (y, x, bottom, right)
    
    rectangles = set()
    
    # Scan for all possible top-left corners
    for y in range(height):
        for x in range(width):
            rect = find_complete_rectangle(y, x)
            if rect:
                rectangles.add(rect)
                
    # Debug: Print rectangle coordinates
    for rect in sorted(rectangles):
        print(f"Rectangle found at: top={rect[0]}, left={rect[1]}, bottom={rect[2]}, right={rect[3]}")
    
    return len(rectangles)

# Input grid (same as before)
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