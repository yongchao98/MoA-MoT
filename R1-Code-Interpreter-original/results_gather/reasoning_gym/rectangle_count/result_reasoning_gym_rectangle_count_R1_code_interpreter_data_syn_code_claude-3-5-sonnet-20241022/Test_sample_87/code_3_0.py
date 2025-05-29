def find_rectangles(grid_str):
    # Convert the grid string into a list of lines
    grid = grid_str.strip().split('\n')
    height = len(grid)
    width = len(grid[0])
    
    def is_valid_rectangle(y1, x1, y2, x2):
        # Check if all corners are border characters
        if not (grid[y1][x1] in '#█' and grid[y1][x2] in '#█' and 
                grid[y2][x1] in '#█' and grid[y2][x2] in '#█'):
            return False
            
        # Check top and bottom edges
        for x in range(x1, x2 + 1):
            if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
                return False
                
        # Check left and right edges
        for y in range(y1, y2 + 1):
            if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
                return False
                
        # Check if it's a minimal rectangle
        inside_chars = set()
        for y in range(y1 + 1, y2):
            for x in range(x1 + 1, x2):
                if grid[y][x] in '#█':
                    inside_chars.add((y, x))
        
        # If there are inside border characters, they must all be '█'
        return all(grid[y][x] == '█' for y, x in inside_chars)
    
    count = 0
    # Find all possible rectangles
    for y1 in range(height-1):
        for y2 in range(y1 + 1, height):
            for x1 in range(width-1):
                for x2 in range(x1 + 1, width):
                    if is_valid_rectangle(y1, x1, y2, x2):
                        count += 1
    
    return count

# Define the test grid
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