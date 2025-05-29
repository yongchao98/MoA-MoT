def find_rectangles(grid):
    # Convert grid to 2D array
    grid = [list(row) for row in grid.splitlines()]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_corner(i, j):
        return grid[i][j] in ['#', '█']
    
    def trace_rectangle(top, left):
        # Find bottom
        bottom = top
        while bottom < height and grid[bottom][left] in ['#', '█']:
            bottom += 1
        bottom -= 1
        
        # Find right
        right = left
        while right < width and grid[top][right] in ['#', '█']:
            right += 1
        right -= 1
        
        # Verify rectangle
        if not all(grid[i][left] in ['#', '█'] for i in range(top, bottom + 1)):
            return None
        if not all(grid[i][right] in ['#', '█'] for i in range(top, bottom + 1)):
            return None
        if not all(grid[top][j] in ['#', '█'] for j in range(left, right + 1)):
            return None
        if not all(grid[bottom][j] in ['#', '█'] for j in range(left, right + 1)):
            return None
            
        return (top, left, bottom, right)
    
    rectangles = set()
    # Scan for top-left corners
    for i in range(height):
        for j in range(width):
            if is_corner(i, j):
                # Check if it's a potential top-left corner
                if (i == 0 or not is_corner(i-1, j)) and (j == 0 or not is_corner(i, j-1)):
                    rect = trace_rectangle(i, j)
                    if rect:
                        rectangles.add(rect)
    
    return len(rectangles)

# Test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                     ###        
                                                                     # #        
            #########################################                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
            #                                       #                # #        
    ########█############                           #                # #        
    #       #           #                           #                # #        
    #       #           #                           #                # #        
    #       #         ##█######################     #                # #        
    #       #         # #                     #     #                # #        
    #       #         # #                     #     #                # #        
    #       #         # #                     #     #                # #        
    #       #         # #                     #     #                # #        
    #       #         # #                     #     #                # #        
    #       #         # #                     #     #                # #        
    #       ##########█#█#####################█######                # #        
    #                 # #                     #                      # #        
    #                 # #                     #                      # #        
    #                 # #                     #                      # #        
    #                 # #                     #                      # #        
    #                 ##█######################                      # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #                   #                                            # #        
    #####################                                            # #        
                                                                     ###        
                                                                                """

print(find_rectangles(grid))