def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    rectangles = 0
    
    # Function to check if a point is part of rectangle border
    def is_border(x, y):
        return 0 <= x < width and 0 <= y < height and grid[y][x] == '#'
    
    # Function to verify if coordinates form a valid rectangle
    def is_valid_rectangle(x1, y1, x2, y2):
        # Check if too small
        if x2 <= x1 or y2 <= y1:
            return False
            
        # Check all border points
        for x in range(x1, x2 + 1):
            if not is_border(x, y1) or not is_border(x, y2):
                return False
        for y in range(y1, y2 + 1):
            if not is_border(x1, y) or not is_border(x2, y):
                return False
        return True
    
    # Find all '#' characters as potential top-left corners
    for y1 in range(height):
        for x1 in range(width):
            if not is_border(x1, y1):
                continue
                
            # For each potential top-left corner, look for bottom-right corner
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if is_border(x2, y2) and is_valid_rectangle(x1, y1, x2, y2):
                        rectangles += 1
                        break  # Once we find a valid rectangle, move to next top-left corner
                if rectangles > rectangles - 1:  # If we found a rectangle
                    break
    
    print(rectangles)

# Create the grid string
grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                      #######                                   
                                      #     #                                   
                                      #     #                                   
                                      #     #                                   
                                      #     #                                   
                                      #######                ###                
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
                                                                                
                                                                                
                                                                                
                                                                                
                                                   ###################          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   #                 #          
                                                   ###################          
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                     #######    
                                                                     #     #    
                                                                     #     #    
                                                                     #     #    
                                                                     #     #    
                                                                     #     #    
                                                                     #######    
"""

find_rectangles(grid)