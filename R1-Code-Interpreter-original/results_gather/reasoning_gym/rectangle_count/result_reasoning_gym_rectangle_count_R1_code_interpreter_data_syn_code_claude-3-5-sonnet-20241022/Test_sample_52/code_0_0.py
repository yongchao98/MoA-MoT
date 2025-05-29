def find_rectangles(grid):
    if not grid:
        return 0
    
    # Convert grid to list of strings for easier processing
    grid = [list(row) for row in grid.strip().split('\n')]
    height = len(grid)
    width = len(grid[0])
    
    def is_corner(i, j):
        # Check if position is a corner (has perpendicular edges)
        if grid[i][j] not in ['#', '█']:
            return False
        
        # Check different corner configurations
        horizontal = ['.', '#', '█']
        vertical = ['.', '#', '█']
        
        # Top-left corner
        if (j == 0 or grid[i][j-1] == ' ') and (i == 0 or grid[i-1][j] == ' '):
            return True
        # Top-right corner
        if (j == width-1 or grid[i][j+1] == ' ') and (i == 0 or grid[i-1][j] == ' '):
            return True
        # Bottom-left corner
        if (j == 0 or grid[i][j-1] == ' ') and (i == height-1 or grid[i+1][j] == ' '):
            return True
        # Bottom-right corner
        if (j == width-1 or grid[i][j+1] == ' ') and (i == height-1 or grid[i+1][j] == ' '):
            return True
            
        # Check for internal corners
        left = j > 0 and grid[i][j-1] in ['#', '█']
        right = j < width-1 and grid[i][j+1] in ['#', '█']
        up = i > 0 and grid[i-1][j] in ['#', '█']
        down = i < height-1 and grid[i+1][j] in ['#', '█']
        
        return (left and up and not right and not down) or \
               (left and down and not right and not up) or \
               (right and up and not left and not down) or \
               (right and down and not left and not up)
    
    # Find all corners
    corners = []
    for i in range(height):
        for j in range(width):
            if is_corner(i, j):
                corners.append((i, j))
    
    # Count rectangles by matching corners
    rectangles = set()
    for i1, j1 in corners:
        for i2, j2 in corners:
            if i2 > i1 and j2 > j1:  # potential opposite corners
                # Check if other two corners exist
                if (i1, j2) in corners and (i2, j1) in corners:
                    # Verify edges
                    valid = True
                    # Check horizontal edges
                    for j in range(j1+1, j2):
                        if grid[i1][j] not in ['#', '█'] or grid[i2][j] not in ['#', '█']:
                            valid = False
                            break
                    # Check vertical edges
                    for i in range(i1+1, i2):
                        if grid[i][j1] not in ['#', '█'] or grid[i][j2] not in ['#', '█']:
                            valid = False
                            break
                    if valid:
                        rectangles.add((i1, j1, i2, j2))
    
    print(len(rectangles))

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

find_rectangles(grid)