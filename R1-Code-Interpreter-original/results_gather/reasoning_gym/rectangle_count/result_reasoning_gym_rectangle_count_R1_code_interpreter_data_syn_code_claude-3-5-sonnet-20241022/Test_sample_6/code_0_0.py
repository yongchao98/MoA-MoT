def find_rectangles(grid):
    if not grid:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    rectangles = set()
    
    def is_valid_corner(char):
        return char in ['#', '█']
    
    def check_rectangle(top, left, bottom, right):
        # Check if this forms a valid rectangle
        if bottom >= rows or right >= cols:
            return False
            
        # Check horizontal borders
        for x in range(left, right + 1):
            if not is_valid_corner(grid[top][x]) or not is_valid_corner(grid[bottom][x]):
                return False
                
        # Check vertical borders
        for y in range(top, bottom + 1):
            if not is_valid_corner(grid[y][left]) or not is_valid_corner(grid[y][right]):
                return False
                
        return True
    
    # Convert input string to 2D array
    grid = [list(row) for row in grid.split('\n')]
    
    # Find all possible rectangles
    for top in range(rows):
        for left in range(cols):
            if is_valid_corner(grid[top][left]):
                # Try all possible bottom-right combinations
                for bottom in range(top, rows):
                    for right in range(left, cols):
                        if check_rectangle(top, left, bottom, right):
                            rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

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