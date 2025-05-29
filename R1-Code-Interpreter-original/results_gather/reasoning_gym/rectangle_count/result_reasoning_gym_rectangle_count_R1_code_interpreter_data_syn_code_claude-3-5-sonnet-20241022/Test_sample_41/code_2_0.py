def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    
    def trace_rectangle(y, x):
        # Find the complete rectangle starting from top-left corner
        if grid[y][x] not in '#█':
            return None
            
        # Find width
        right = x
        while right < width and grid[y][right] in '#█':
            right += 1
        right -= 1
        
        # Find height
        bottom = y
        while bottom < height and grid[bottom][x] in '#█':
            bottom += 1
        bottom -= 1
        
        # Verify it's a complete rectangle
        # Check all borders
        for i in range(x, right + 1):
            if grid[y][i] not in '#█' or grid[bottom][i] not in '#█':
                return None
        for i in range(y, bottom + 1):
            if grid[i][x] not in '#█' or grid[i][right] not in '#█':
                return None
                
        # Check interior (must have some empty space)
        has_space = False
        for i in range(y + 1, bottom):
            for j in range(x + 1, right):
                if grid[i][j] == ' ':
                    has_space = True
                    break
            if has_space:
                break
        
        if not has_space:
            return None
            
        return (x, y, right, bottom)

    rectangles = set()
    # Scan for top-left corners of rectangles
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Check if this could be a top-left corner
                if (y == 0 or grid[y-1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█'):
                    rect = trace_rectangle(y, x)
                    if rect:
                        rectangles.add(rect)

    return len(rectangles)

grid = """                                                                                
                                                                                
                                                                                
                                #############################                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                         #######█#######################    #                   
                         #      #                      #    #                   
                         #      #                      #    #                   
                         #      #                      #    #                   
                         #######█#######################    #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #                           #                   
                                #############################                   
                                                                                
                                                                                
               ###################################                              
               #                                 #                              
               #                                 #                              
               #                                 #                              
               #                                 #                              
               #                                 #            ###############   
               #                                 #            #             #   
               ###################################            #             #   
                                                              #             #   
                                                              #             #   
                                                              #             #   
                                                              #             #   
                                                              #             #   
                                                              #             #   
                                                           ###█########     #   
                                                           #  ########█######   
                                                           ############         
                                                                                
                                                                                """

print(find_rectangles(grid))