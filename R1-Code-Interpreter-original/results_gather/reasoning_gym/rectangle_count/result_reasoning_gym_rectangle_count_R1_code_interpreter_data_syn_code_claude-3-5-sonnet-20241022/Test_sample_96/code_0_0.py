def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    grid = [line for line in grid.split('\n') if line.strip()]
    
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    rectangles = 0
    
    # Function to check if a point is part of rectangle border
    def is_border(x, y):
        return 0 <= x < width and 0 <= y < height and grid[y][x] in '#█'
    
    # Function to verify if coordinates form a valid rectangle
    def is_valid_rectangle(x1, y1, x2, y2):
        # Check corners
        if not all(is_border(x, y) for x, y in [(x1, y1), (x2, y1), (x1, y2), (x2, y2)]):
            return False
            
        # Check horizontal borders
        for x in range(x1 + 1, x2):
            if not (is_border(x, y1) and is_border(x, y2)):
                return False
                
        # Check vertical borders
        for y in range(y1 + 1, y2):
            if not (is_border(x1, y) and is_border(x2, y)):
                return False
                
        return True
    
    # Find all '#' or '█' characters as potential corners
    corners = [(x, y) for y in range(height) for x in range(width) if is_border(x, y)]
    
    # Check all possible rectangle combinations
    found_rectangles = set()
    for i, (x1, y1) in enumerate(corners):
        for x2, y2 in corners[i+1:]:
            if x2 > x1 and y2 > y1:  # Ensure second point is below and right of first point
                if is_valid_rectangle(x1, y1, x2, y2):
                    found_rectangles.add((x1, y1, x2, y2))
    
    print(len(found_rectangles))

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