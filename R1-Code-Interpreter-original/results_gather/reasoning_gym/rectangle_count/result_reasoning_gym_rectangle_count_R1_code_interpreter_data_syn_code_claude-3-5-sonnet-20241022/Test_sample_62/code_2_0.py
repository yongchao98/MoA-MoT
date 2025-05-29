def find_rectangles(grid):
    # Clean and convert grid to 2D array
    grid = [list(line) for line in grid.strip().splitlines()]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_rectangle_corner(i, j):
        # Check if position is valid
        if i < 0 or i >= height or j < 0 or j >= width:
            return False
        # Check if it's a corner character
        if grid[i][j] not in ['#', '█']:
            return False
            
        # Count adjacent edges
        horizontal = 0
        vertical = 0
        
        # Check horizontal neighbors
        if j > 0 and grid[i][j-1] in ['#', '█']:
            horizontal += 1
        if j < width-1 and grid[i][j+1] in ['#', '█']:
            horizontal += 1
            
        # Check vertical neighbors
        if i > 0 and grid[i-1][j] in ['#', '█']:
            vertical += 1
        if i < height-1 and grid[i+1][j] in ['#', '█']:
            vertical += 1
            
        # It's a corner if we have exactly one horizontal and one vertical connection
        return horizontal >= 1 and vertical >= 1

    def trace_rectangle(start_i, start_j):
        # Find the dimensions of the rectangle
        right_j = start_j
        bottom_i = start_i
        
        # Find right edge
        while right_j < width and grid[start_i][right_j] in ['#', '█']:
            right_j += 1
        right_j -= 1
        
        # Find bottom edge
        while bottom_i < height and grid[bottom_i][start_j] in ['#', '█']:
            bottom_i += 1
        bottom_i -= 1
        
        # Verify rectangle
        for i in range(start_i, bottom_i + 1):
            if grid[i][right_j] not in ['#', '█']:
                return None
        for j in range(start_j, right_j + 1):
            if grid[bottom_i][j] not in ['#', '█']:
                return None
        
        return (start_i, start_j, bottom_i, right_j)

    rectangles = set()
    visited = set()
    
    # Scan for corners
    for i in range(height):
        for j in range(width):
            if grid[i][j] in ['#', '█'] and (i, j) not in visited:
                if is_rectangle_corner(i, j):
                    rect = trace_rectangle(i, j)
                    if rect:
                        rectangles.add(rect)
                        # Mark rectangle points as visited
                        for y in range(rect[0], rect[2] + 1):
                            for x in range(rect[1], rect[3] + 1):
                                if grid[y][x] in ['#', '█']:
                                    visited.add((y, x))
    
    return len(rectangles)

# Test grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                           #############################################        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                                           #        
                           #                               ########    #        
                           #                               #      #    #        
                       ####█#####                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          #      #    #        
                       #   #    #                          ########    #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
                       #   #    #                                      #        
         ##############█###█####█##########################         ####        
         #             ####█#####                         #         # ##        
         #                 #                              #         ####        
         #                 #                           ###█##########  #        
         #                 ############################█##█#########█###        
         # #########################################   #  #         #           
         # #                                       #   #  #         #           
         # #                                       #   #  #         #           
         # #                                       #   █###         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #                                       #   █ ##         #           
         # #########################################   █###         #           
         #                                             ###█##########           
         #                                                #                     
         #                                                #                     
         #                                                #                     
         #                                                #                     
         #                             ###################█##########           
         ##############################█###################         #           
                                       #                            #           
                                       ##############################           """

print(find_rectangles(grid))