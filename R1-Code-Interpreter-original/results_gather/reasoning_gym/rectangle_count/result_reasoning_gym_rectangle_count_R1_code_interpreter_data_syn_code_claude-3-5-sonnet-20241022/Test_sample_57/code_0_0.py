def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.splitlines() if line.strip()]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    count = 0
    
    # Function to check if a rectangle is valid
    def is_valid_rectangle(top, left, bottom, right):
        # Check corners
        if grid[top][left] not in '#█' or grid[top][right] not in '#█' or \
           grid[bottom][left] not in '#█' or grid[bottom][right] not in '#█':
            return False
            
        # Check horizontal lines
        for x in range(left + 1, right):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
                
        # Check vertical lines
        for y in range(top + 1, bottom):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
                
        return True
    
    # Find all possible rectangles
    for top in range(height):
        for bottom in range(top + 1, height):
            for left in range(width):
                for right in range(left + 1, width):
                    if is_valid_rectangle(top, left, bottom, right):
                        count += 1
    
    print(count)

# Create the grid string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
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
                                                                             # #
                                                                             # #
                                                                             # #
                                                                             # #
                                                                             # #
                                                                             # #
                                                                             ###
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
        ######################                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        #                    #                                                  
        ######################                                                  
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

find_rectangles(grid)