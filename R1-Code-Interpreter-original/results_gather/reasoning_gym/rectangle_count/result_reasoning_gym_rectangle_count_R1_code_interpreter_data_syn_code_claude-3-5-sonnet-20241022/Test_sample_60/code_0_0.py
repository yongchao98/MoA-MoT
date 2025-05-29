def find_rectangles(grid):
    if not grid:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    
    def is_border(char):
        return char in '#█'
    
    def verify_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        
        return True
    
    rectangles = set()
    
    # Find all potential top-left corners
    for y in range(height):
        for x in range(width):
            if is_border(grid[y][x]):
                # Look for bottom-right corners
                for bottom in range(y + 1, height):
                    for right in range(x + 1, width):
                        if is_border(grid[bottom][right]):
                            if verify_rectangle(y, x, bottom, right):
                                rectangles.add((y, x, bottom, right))
    
    return len(rectangles)

# Convert the ASCII art into a list of strings
grid_str = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
         ##############################################                         
         #                                            #                         
         #                                            #                         
         #                                            #                         
         #                                            #                         
         #                                            #                         
         #                                            #                         
         ##############################################                         
                     ################################################           
                     #                                              #           
                     #                                              #           
                     #                                              #           
                     #                                              #           
            #########█##############################################█######     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #     
            #        #                                              #     #  ###
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #        #                                              #     #  # #
            #      ##█##############################################█###  #  # #
            #      # ################################################  #  #  # #
            #      #                                                   #  #  # #
            #      #                                                   #  #  # #
            #      #####################################################  #  # #
            #                                                             #  # #
            #                                                             #  # #
            #                                                             #  # #
            #                                                             #  # #
            #                                                             #  # #
            #                                                             #  ###
            #                                                     ########█###  
            #                                                     #       #  #  
            #                                                     #       #  #  
            #                                                     #       #  #  
   #########█##############                                       #       #  #  
   #        ##############█#######################################█########  #  
   #                      #                                       ############  
   ########################                                                     
                                                                                
                                                                                
"""

# Convert the string into a list of lines
grid = [list(line) for line in grid_str.splitlines()]

# Remove empty lines at the beginning and end
while grid and not any(grid[0]):
    grid.pop(0)
while grid and not any(grid[-1]):
    grid.pop()

print(find_rectangles(grid))