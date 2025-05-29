def find_rectangles(grid_str):
    # Convert the grid string into a list of lines
    grid = [line for line in grid_str.strip().split('\n')]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_border(char):
        return char in '#█'
    
    def check_rectangle(top, left, bottom, right):
        # Check if all edges form a valid rectangle
        # Top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        # Check if it's a minimal rectangle (no other rectangles inside)
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if is_border(grid[y][x]) and grid[y][x] != '█':
                    return False
        return True

    rectangles = set()
    # Find all possible rectangles
    for top in range(height):
        for bottom in range(top + 1, height):
            for left in range(width):
                for right in range(left + 1, width):
                    if check_rectangle(top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

# Define the grid as a multi-line string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                       ######                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                               ##############      ####█####█###############    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               #            #      #   #    #              #    
                               ##############      #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   #   #    #              #    
                                                   ####█####█###############    
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       #    #                   
                                                       ######                   
                                                                                
                       ##################################################       
                       #                                                #       
                       #                                                #       
                       #                                                #       
                       #                                                #       
                       #                                                #       
                       #  #################################             #       
                       #  #                               #             #       
                       #  #################################             #       
                       #                                                #       
                       ##################################################       
                                                                                """

print(find_rectangles(grid))