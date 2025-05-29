def find_rectangles(grid):
    # Convert grid to 2D array
    grid = [list(row) for row in grid.strip('\n').split('\n')]
    height = len(grid)
    width = len(grid[0])
    
    def is_border_char(c):
        return c in '#█'
    
    def validate_rectangle(top, left, bottom, right):
        # Check if this forms a valid rectangle
        if top >= bottom or left >= right:
            return False
            
        # Check all borders
        for x in range(left, right + 1):
            if not is_border_char(grid[top][x]) or not is_border_char(grid[bottom][x]):
                return False
        for y in range(top, bottom + 1):
            if not is_border_char(grid[y][left]) or not is_border_char(grid[y][right]):
                return False
                
        return True
    
    rectangles = set()
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if not is_border_char(grid[top][left]):
                continue
                
            for bottom in range(top + 1, height):
                if not is_border_char(grid[bottom][left]):
                    continue
                    
                for right in range(left + 1, width):
                    if not is_border_char(grid[bottom][right]) or not is_border_char(grid[top][right]):
                        continue
                        
                    if validate_rectangle(top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
    # Count overlapping rectangles
    overlaps = set()
    for r1 in rectangles:
        for r2 in rectangles:
            if r1 >= r2:
                continue
            # Check if rectangles overlap
            if not (r1[2] < r2[0] or r2[2] < r1[0] or  # no vertical overlap
                   r1[3] < r2[1] or r2[3] < r1[1]):    # no horizontal overlap
                overlaps.add((r1, r2))
    
    # Each '█' character represents one overlap between two rectangles
    # Total rectangles = number of distinct rectangles found
    return len(rectangles)

# The grid as a string
grid = """                                                                                
                                      ####                                      
          ############################█##█#############################         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          #                           #  #                            #         
          ############################█##█#############################         
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      #  #                                      
                                      ####                                      
                                                                                
                                                                                
                                                ##################              
                                                #                #              
                                                #                #              
                                                #                #              
                                                #                #              
                                                #                #              
                    ############################██████████████████###           
                    #                                               #           
                    #                                               #           
                    #                                               #           
                    #################################################           
                                                                                """

print(find_rectangles(grid))