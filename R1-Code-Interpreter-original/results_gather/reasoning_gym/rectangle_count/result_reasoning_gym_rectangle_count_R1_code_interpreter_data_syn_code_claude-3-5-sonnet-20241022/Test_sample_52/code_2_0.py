def count_rectangles(grid):
    # Convert grid to list of lines and clean up
    lines = [line.rstrip() for line in grid.strip().split('\n')]
    
    # Find dimensions
    height = len(lines)
    width = max(len(line) for line in lines)
    
    # Normalize grid width and create 2D array
    grid = [[' ' if j >= len(lines[i]) else lines[i][j] 
             for j in range(width)] for i in range(height)]

    def trace_rectangle(i, j):
        # Find right edge
        right = j
        while right < width and grid[i][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = i
        while bottom < height and grid[bottom][j] in '#█':
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        if right <= j or bottom <= i:
            return False
            
        # Check all corners
        if grid[i][j] not in '#█' or grid[i][right] not in '#█' or \
           grid[bottom][j] not in '#█' or grid[bottom][right] not in '#█':
            return False
            
        # Check top and bottom edges
        for x in range(j, right + 1):
            if grid[i][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
                
        # Check left and right edges
        for y in range(i, bottom + 1):
            if grid[y][j] not in '#█' or grid[y][right] not in '#█':
                return False
                
        # Check that we don't have any other '#' or '█' touching the edges
        if j > 0:
            for y in range(i, bottom + 1):
                if grid[y][j-1] in '#█':
                    return False
        if right < width - 1:
            for y in range(i, bottom + 1):
                if grid[y][right+1] in '#█':
                    return False
        if i > 0:
            for x in range(j, right + 1):
                if grid[i-1][x] in '#█':
                    return False
        if bottom < height - 1:
            for x in range(j, right + 1):
                if grid[bottom+1][x] in '#█':
                    return False
                    
        return True

    # Count rectangles by finding top-left corners
    count = 0
    rectangles_found = set()
    
    for i in range(height):
        for j in range(width):
            # Look for potential top-left corners
            if grid[i][j] in '#█':
                # Check if it's a top-left corner
                if (j == 0 or grid[i][j-1] not in '#█') and \
                   (i == 0 or grid[i-1][j] not in '#█'):
                    if trace_rectangle(i, j):
                        count += 1

    print(count)

# Your ASCII grid
grid = """    #################################################################           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                                               #           
    #                                              #################█########   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                              #                #       #   
    #                                   ###########████#############█########   
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             #           
    #                                   #             #             ########    
    ####################################█#############█###############     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        #             #              #     #    
                                        ###############              #     #    
                                                                     #######    
                                                                      ########  
                                                                      #      #  
                                                                      ########  """

count_rectangles(grid)