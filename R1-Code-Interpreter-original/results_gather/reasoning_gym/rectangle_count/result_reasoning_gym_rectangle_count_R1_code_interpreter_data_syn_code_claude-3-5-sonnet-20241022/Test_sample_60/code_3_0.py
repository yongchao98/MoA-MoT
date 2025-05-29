def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_corner(y, x):
        # Check if position is a valid corner character
        if grid[y][x] not in '#█':
            return False
            
        # For top-left corner:
        # No border character above and no border character to the left
        above = grid[y-1][x] if y > 0 else ' '
        left = grid[y][x-1] if x > 0 else ' '
        return above not in '#█' and left not in '#█'
    
    def trace_rectangle(start_y, start_x):
        # Find right end of top edge
        right = start_x + 1
        while right < width and grid[start_y][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom end of left edge
        bottom = start_y + 1
        while bottom < height and grid[bottom][start_x] in '#█':
            bottom += 1
        bottom -= 1
        
        # Verify right edge
        for y in range(start_y, bottom + 1):
            if grid[y][right] not in '#█':
                return False
                
        # Verify bottom edge
        for x in range(start_x, right + 1):
            if grid[bottom][x] not in '#█':
                return False
                
        return (bottom, right)
    
    rectangles = set()
    
    # Find all top-left corners and trace rectangles
    for y in range(height):
        for x in range(width):
            if is_corner(y, x):
                result = trace_rectangle(y, x)
                if result:
                    bottom, right = result
                    # Store rectangle coordinates
                    rectangles.add((y, x, bottom, right))
    
    return len(rectangles)

# Input grid
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

# Clean up the grid and remove empty lines
grid = [list(line) for line in grid_str.splitlines() if line.strip()]

print(find_rectangles(grid))