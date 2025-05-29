def find_rectangles(grid):
    # Convert grid to 2D array and clean up
    lines = [line for line in grid.splitlines() if line.strip()]
    grid = [list(line) for line in lines]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0

    def is_top_left_corner(i, j):
        if grid[i][j] not in ['#', '█']:
            return False
            
        # Check if it's a top-left corner
        top = i == 0 or grid[i-1][j] not in ['#', '█']
        left = j == 0 or grid[i][j-1] not in ['#', '█']
        right = j < width-1 and grid[i][j+1] in ['#', '█']
        bottom = i < height-1 and grid[i+1][j] in ['#', '█']
        
        return top and left and right and bottom

    def validate_rectangle(i, j):
        # Find right edge
        right = j
        while right < width and grid[i][right] in ['#', '█']:
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = i
        while bottom < height and grid[bottom][j] in ['#', '█']:
            bottom += 1
        bottom -= 1
        
        # Verify all edges
        for x in range(j, right + 1):
            if grid[i][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
                return None
        for y in range(i, bottom + 1):
            if grid[y][j] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
                return None
                
        return (i, j, bottom, right)

    rectangles = set()
    
    # Find all rectangles starting from top-left corners
    for i in range(height):
        for j in range(width):
            if is_top_left_corner(i, j):
                rect = validate_rectangle(i, j)
                if rect:
                    rectangles.add(rect)
                    
                    # If this corner is a '█', look for potential overlapping rectangle
                    if grid[i][j] == '█':
                        # Check for additional rectangle starting from adjacent positions
                        if i > 0 and grid[i-1][j] in ['#', '█']:
                            upper_rect = validate_rectangle(i-1, j)
                            if upper_rect:
                                rectangles.add(upper_rect)
                        if j > 0 and grid[i][j-1] in ['#', '█']:
                            left_rect = validate_rectangle(i, j-1)
                            if left_rect:
                                rectangles.add(left_rect)

    return len(rectangles)

# Test grid
test_grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
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

print(find_rectangles(test_grid))