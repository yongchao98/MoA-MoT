def find_rectangles(grid):
    if not grid:
        return 0
    
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line.strip()]
    
    height = len(lines)
    width = len(lines[0])
    
    # Convert grid to 2D array for easier processing
    grid = [[c for c in line] for line in lines]
    
    def is_border(c):
        return c in '#█'
    
    def check_rectangle(start_y, start_x):
        if not is_border(grid[start_y][start_x]):
            return False
            
        # Find right edge
        right = start_x
        while right < width and is_border(grid[start_y][right]):
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = start_y
        while bottom < height and is_border(grid[bottom][start_x]):
            bottom += 1
        bottom -= 1
        
        # Verify rectangle
        for y in range(start_y, bottom + 1):
            if not is_border(grid[y][start_x]) or not is_border(grid[y][right]):
                return False
        for x in range(start_x, right + 1):
            if not is_border(grid[start_y][x]) or not is_border(grid[bottom][x]):
                return False
                
        # Mark as visited
        for y in range(start_y, bottom + 1):
            for x in range(start_x, right + 1):
                if grid[y][x] == '#':
                    grid[y][x] = '.'
                elif grid[y][x] == '█':
                    grid[y][x] = '#'
        
        return True
    
    count = 0
    # First pass: count regular rectangles and mark overlapping ones
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '#':
                if check_rectangle(y, x):
                    count += 1
    
    # Second pass: count remaining rectangles (from overlapping parts)
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '#':
                if check_rectangle(y, x):
                    count += 1
    
    return count

# Read the grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                               #########        
                                                               #       #        
                                                               #       #        
                                                               #       #        
                                                              #█#      #        
                                                              ###      #        
                                                              ###      #        
                                                              ###      #        
                                                              ###      #        
                                                              ###      #        
                                                              ###      #        
                                                              ###      #        
                                                              ###      #        
                    ##########################################███#     #        
                    #                                         ####     #        
                    #                                         ####     #        
                    #                                         ####     #        
                    #                                         ####     #        
                    #                                         ####     #        
                    #                                         ####     #        
                    #                                         ####     #        
                  ##█######################                   ####     #        
                  # #                     #                   ####     #        
                  # # ##################  #                   ####     #        
                  # # #                #  #                   ####     #        
                  # # #                #  #                   ####     #        
                  # # #                #  #                   ####     #        
       ###########█#█#█################█##█###################████##   #        
       #          # # #                #  #                   #### #   #        
       #          # # #                #  #                   #### #   #        
       #          # # #                #  #                   #### #   #        
       #          # # #                #  #         ##########████#█###█########
       #          # # #                #  #         #         #### #   #       #
       #          # # #                #  #         #         ████##   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         █#####   #       #
       #          # # #                #  #         #         ██####   #       #
       #          # # #                #  #         #         ## ###   #       #
       #          # # #                #  #         #         ## ###   #       #
       #          # # #                #  #         #         ## ###   #       #
       #          # # #                #  #         #         ## ###   #       #
       #          # # #                #  #         #         ## ###   #       #
       #          # # ##################  #     ####█#####    ## ###   #       #
       #          # #                     #     #   #    #    ## ###   #       #
       #          # #                     #     #   #    #    ## ###   #       #
       #          # #                     #     ####█#####    ## ###   #       #
       #          # #                     #         #         ## ###   #       #
       #          # #                     #         #         ###███####       #
       #          # # ####################█#########█##       #  ###           #
       #          # # #                   #         # #       #  ###           #
       #          # # #                   #         # #       #  ###           #
       #          # # #                   #         # #       #  ###           #
       #          # # #                   #         # #       #  ###           #
       #          # # #                   #         # #       #  ###           #
       #          # # #                   #         # #       #  ###           #
       #          # # #                   #         ##█#######█##███############
       #          # # #                   #           #       #  ###            
       #          # # ####################█############       #  ###            
       #          # #                     #                   #  ###            
       #          # #                     #                   #  ###            
       #          # #                     #                   #  ###            
       #          # #                     #                   ###█##            
       ###########█#█#####################█######################█##            
                  ##█######################                      #              
                    #                                            #              
                    ##############################################              
                                                                                
                                                                                
"""

print(find_rectangles(grid))