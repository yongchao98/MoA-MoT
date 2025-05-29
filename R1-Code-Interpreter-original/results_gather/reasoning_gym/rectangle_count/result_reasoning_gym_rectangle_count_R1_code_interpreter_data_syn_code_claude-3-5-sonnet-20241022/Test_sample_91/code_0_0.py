def count_rectangles(grid):
    # Convert the grid into a list of strings for easier processing
    grid = [line for line in grid.splitlines()]
    
    height = len(grid)
    if height == 0:
        return 0
    width = len(grid[0])
    
    # Function to check if a point is part of a rectangle border
    def is_border(x, y):
        return 0 <= x < width and 0 <= y < height and grid[y][x] in '#█'
    
    # Function to validate if points form a valid rectangle
    def is_valid_rectangle(x1, y1, x2, y2):
        # Check all corners
        if not (is_border(x1, y1) and is_border(x2, y1) and 
                is_border(x1, y2) and is_border(x2, y2)):
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
    
    rectangles = set()
    # Scan for rectangles
    for y1 in range(height):
        for x1 in range(width):
            if is_border(x1, y1):
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if is_valid_rectangle(x1, y1, x2, y2):
                            rectangles.add((x1, y1, x2, y2))
    
    print(len(rectangles))

# Create the grid string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                  ########                                      
                                  #      #                                      
                                  #      #                                      
                                  #      #                                      
                                  #      #                                      
                                  #      #                                      
                                  #      #                                      
                                  #      #                                      
                                  #      #                                      
                                  #      #                                      
                                  ########                                      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                """

count_rectangles(grid)