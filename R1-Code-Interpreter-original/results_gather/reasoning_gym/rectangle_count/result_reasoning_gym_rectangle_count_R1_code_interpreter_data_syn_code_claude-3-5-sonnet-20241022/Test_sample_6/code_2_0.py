def find_rectangles(grid):
    # Split the grid into lines
    grid = [list(line) for line in grid.split('\n')]
    
    if not grid:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    rectangles = set()
    
    def is_valid_corner(char):
        return char in ['#', '█']
    
    def is_empty_inside(top, left, bottom, right):
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if grid[y][x] not in [' ', '#', '█']:
                    return False
        return True
    
    def check_rectangle(top, left, bottom, right):
        # Check if corners are valid
        if not (is_valid_corner(grid[top][left]) and 
                is_valid_corner(grid[top][right]) and 
                is_valid_corner(grid[bottom][left]) and 
                is_valid_corner(grid[bottom][right])):
            return False
        
        # Check horizontal borders
        for x in range(left + 1, right):
            if not is_valid_corner(grid[top][x]):
                return False
            if not is_valid_corner(grid[bottom][x]):
                return False
        
        # Check vertical borders
        for y in range(top + 1, bottom):
            if not is_valid_corner(grid[y][left]):
                return False
            if not is_valid_corner(grid[y][right]):
                return False
        
        # Check if inside is valid
        return is_empty_inside(top, left, bottom, right)
    
    # Find top-left corners
    for top in range(rows):
        for left in range(cols):
            if is_valid_corner(grid[top][left]):
                # Find matching bottom-right corners
                for bottom in range(top + 1, rows):
                    for right in range(left + 1, cols):
                        if (grid[top][left] == '#' or grid[top][left] == '█') and \
                           check_rectangle(top, left, bottom, right):
                            # Found a valid rectangle
                            rectangles.add((top, left, bottom, right))
    
    # Count distinct rectangles
    count = 0
    processed = set()
    
    for rect in rectangles:
        top, left, bottom, right = rect
        # Check if this rectangle is unique
        if (top, left, bottom, right) not in processed:
            count += 1
            processed.add((top, left, bottom, right))
    
    return count

# Input grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                      ##########################
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      #                        #
                                                      ##########################
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                       ###################      
                                                       #                 #      
                                                       ###################      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                         ####################                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #    ##############█#################  
                                         #    #             #                #  
                                         #    #             #                #  
                                         #    ##############█#################  
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
                                         #                  #                   
          ###############################█##################█###################
          #                              ####################                  #
          #                                                                    #
          #                                                                    #
          #                                                                    #
          #                                                                    #
          ######################################██████████████##################
                                                #            #                  
                                                #            #                  
                                                #            #                  
                                                #            #                  
                                                #            #                  
                                                #            #                  
                                                #            #               ###
                                                #            #               # #
                                                #            #               # #
                                                #            #               ###
                                                ##############                  """

print(find_rectangles(grid))