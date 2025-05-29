def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    height = len(grid)
    width = len(grid[0])
    rectangles = 0

    # Function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'

    # Function to validate if coordinates form a valid rectangle
    def is_valid_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True

    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if is_border(grid[top][left]):
                # Try all possible bottom-right corners
                for bottom in range(top, height):
                    for right in range(left, width):
                        if (is_border(grid[top][left]) and 
                            is_border(grid[top][right]) and 
                            is_border(grid[bottom][left]) and 
                            is_border(grid[bottom][right])):
                            if is_valid_rectangle(top, left, bottom, right):
                                # Check if it's a proper rectangle (not just a line)
                                if bottom > top and right > left:
                                    rectangles += 1

    print(rectangles)

# Your ASCII grid
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                              ############      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                     #########################################█####      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #                                        #   #      #      
                     #########################################█####      #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
    ###########################################################          #      
    #                                                        ##          #      
    #                                                        ##          #      
    #                                                        ##          #      
    #                                                        ##          #      
    ###########################################################          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              #          #      
                                                              ############      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

find_rectangles(grid)