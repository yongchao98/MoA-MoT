def find_rectangles(grid_str):
    # Convert the input string into a list of lines
    grid = [line for line in grid_str.split('\n') if line]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_border(char):
        return char in '#█'
    
    def find_complete_rectangle(start_y, start_x):
        # Find the bottom-right corner of the rectangle
        # First, find the right edge
        right_x = start_x
        while right_x < width and is_border(grid[start_y][right_x]):
            right_x += 1
        right_x -= 1
        
        # Then find the bottom edge
        bottom_y = start_y
        while bottom_y < height and is_border(grid[bottom_y][start_x]):
            bottom_y += 1
        bottom_y -= 1
        
        # Validate the rectangle
        # Check if all corners exist
        if not (is_border(grid[start_y][start_x]) and 
                is_border(grid[start_y][right_x]) and
                is_border(grid[bottom_y][start_x]) and
                is_border(grid[bottom_y][right_x])):
            return None
        
        # Check all edges
        for x in range(start_x, right_x + 1):
            if not (is_border(grid[start_y][x]) and is_border(grid[bottom_y][x])):
                return None
        for y in range(start_y, bottom_y + 1):
            if not (is_border(grid[y][start_x]) and is_border(grid[y][right_x])):
                return None
                
        return (start_y, start_x, bottom_y, right_x)
    
    rectangles = set()
    # Find all rectangles
    for y in range(height):
        for x in range(width):
            if is_border(grid[y][x]):
                # Check if this could be the top-left corner of a rectangle
                if (x == 0 or not is_border(grid[y][x-1])) and (y == 0 or not is_border(grid[y-1][x])):
                    rect = find_complete_rectangle(y, x)
                    if rect:
                        rectangles.add(rect)
    
    # Count overlapping points
    overlap_count = 0
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                # Count how many rectangles contain this point
                containing_rects = sum(1 for rect in rectangles 
                                     if rect[0] <= y <= rect[2] and rect[1] <= x <= rect[3])
                if containing_rects > 1:
                    overlap_count += 1
    
    return len(rectangles)

grid_str = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                       #############                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                       #           #                                            
                    ###█###########█######################                      
                    #  #############                     #                      
                    #                                    #                      
                    #                                    #                      
                    #                                    #                      
                    #                                    #                      
                    #                     ###############█#########             
                    #                     #              #        #             
                    #                     #              #        #             
                    #                     ###############█#########             
                    #                                    #                      
                    #                                    #                      
                    #                                    #                      
                    #                                    #                      
                    #                                    #                      
               #####█####################################█#####                 
               #    #                                    #    #                 
               #    #                                    #    #                 
               #    #                                    #    #                 
               #    #                                    #    #                 
               #    #                                    #    #                 
               #    #                                    #    #                 
               #    #                                    #    #                 
               #    #                                    #    #                 
               #    ######################################    #                 
               #                                              #                 
               ################################################                 
"""

print(find_rectangles(grid_str))