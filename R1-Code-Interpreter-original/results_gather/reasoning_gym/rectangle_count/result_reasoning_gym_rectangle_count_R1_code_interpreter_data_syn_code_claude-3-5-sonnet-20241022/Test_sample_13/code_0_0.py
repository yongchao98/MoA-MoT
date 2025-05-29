def find_rectangles(grid):
    if not grid:
        return 0
    
    # Convert string grid to 2D array
    grid = [list(row) for row in grid.strip().split('\n')]
    height = len(grid)
    width = len(grid[0])
    
    def is_rectangle(top, left, bottom, right):
        # Check if the given coordinates form a valid rectangle
        if top >= bottom or left >= right:
            return False
            
        # Check horizontal borders
        for x in range(left, right + 1):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
                
        # Check vertical borders
        for y in range(top, bottom + 1):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
                
        return True
    
    rectangles = 0
    # Check each possible rectangle
    for top in range(height):
        for left in range(width):
            if grid[top][left] not in '#█':
                continue
                
            for bottom in range(top, height):
                for right in range(left, width):
                    if grid[bottom][right] not in '#█':
                        continue
                    if is_rectangle(top, left, bottom, right):
                        # Count overlapping rectangles
                        has_overlap = any(grid[y][x] == '█' 
                                        for y in range(top, bottom + 1)
                                        for x in range(left, right + 1))
                        rectangles += 2 if has_overlap else 1
                        
    return rectangles

# Test with the provided ASCII art
ascii_art = """                                                                                
                                                                                
                                      ########################                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      #                      #                  
                                      ########################                  
                                                             ###                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                     ########█#█#######         
                                                     #       # #      #         
        ####################                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         #       # #      #         
        #                  #                         ########█#█#######         
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        #                  #                                 # #                
        ####################                                 # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             # #                
                                                             ###                
                                                                                
                           ###################################################  
                           #                                                 #  
                           #                                                 #  
                           #                                                 #  
                           ###################################################  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

print(find_rectangles(ascii_art))