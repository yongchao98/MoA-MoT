def find_rectangles(grid):
    # Convert grid to list of strings
    lines = grid.split('\n')
    
    # Find first and last non-empty lines to get actual grid boundaries
    start_line = 0
    end_line = len(lines) - 1
    
    while start_line < len(lines) and not lines[start_line].strip():
        start_line += 1
    while end_line >= 0 and not lines[end_line].strip():
        end_line -= 1
        
    # Get actual grid
    grid = lines[start_line:end_line+1]
    
    # Count rectangles by identifying their top-left corners and following their borders
    rectangles = 0
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_top_left_corner(y, x):
        if y < 0 or x < 0 or y >= height or x >= width:
            return False
        current = grid[y][x]
        if current not in '#█':
            return False
        # Check if it's a corner by looking at adjacent positions
        left_empty = x == 0 or grid[y][x-1] == ' '
        top_empty = y == 0 or grid[y-1][x] == ' '
        return left_empty and top_empty
    
    def validate_rectangle(y, x):
        if not is_top_left_corner(y, x):
            return False
            
        # Find right edge
        right = x
        while right < width and grid[y][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = y
        while bottom < height and grid[bottom][x] in '#█':
            bottom += 1
        bottom -= 1
        
        # Validate complete rectangle
        if bottom - y < 1 or right - x < 1:
            return False
            
        # Check all borders
        for i in range(x, right + 1):
            if grid[y][i] not in '#█' or grid[bottom][i] not in '#█':
                return False
        for i in range(y, bottom + 1):
            if grid[i][x] not in '#█' or grid[i][right] not in '#█':
                return False
                
        return True
    
    # Scan for rectangles
    for y in range(height):
        for x in range(width):
            if validate_rectangle(y, x):
                rectangles += 1
                
    print(rectangles)

# Your ASCII grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                              ############      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                     #########################################█####      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #########################################█####      #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
    ###########################################################          #      
    #                                                        ##          #      
    #                                                        ##          #      
    #                                                        ##          #      
    #                                                        ##          #      
    ###########################################################          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              ############      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

find_rectangles(grid)