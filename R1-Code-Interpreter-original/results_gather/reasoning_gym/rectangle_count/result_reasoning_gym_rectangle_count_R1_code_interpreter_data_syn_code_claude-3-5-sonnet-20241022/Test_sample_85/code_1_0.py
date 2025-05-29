def find_rectangles(grid_str):
    # Convert the input string into a list of lines
    grid = [line for line in grid_str.split('\n') if line]
    
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Helper function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'
    
    # Helper function to validate rectangle
    def validate_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True
    
    # Find all potential rectangles
    for y1 in range(height):
        for x1 in range(width):
            if is_border(grid[y1][x1]):
                # Found potential top-left corner
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if is_border(grid[y2][x2]):
                            # Check if this forms a valid rectangle
                            if validate_rectangle(y1, x1, y2, x2):
                                # Verify it's not just a line
                                if y2 > y1 and x2 > x1:
                                    rectangles.append((y1, x1, y2, x2))
    
    # Remove duplicates
    rectangles = list(set(rectangles))
    
    # Count overlapping rectangles
    overlaps = set()
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                count = 0
                for rect in rectangles:
                    top, left, bottom, right = rect
                    if left <= x <= right and top <= y <= bottom:
                        count += 1
                if count > 1:
                    overlaps.add((y, x))
    
    return len(rectangles) + len(overlaps)

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