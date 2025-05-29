def find_rectangles(grid):
    # Convert grid to list of strings for easier processing
    grid = [line for line in grid.split('\n') if line.strip()]
    
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Function to check if a point is part of rectangle border
    def is_border(char):
        return char in '#█'
    
    # Function to validate if points form a valid rectangle
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
    for y1 in range(height):
        for x1 in range(width):
            if is_border(grid[y1][x1]):
                # Found potential top-left corner
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if is_border(grid[y2][x2]):
                            # Check if these points form a valid rectangle
                            if is_valid_rectangle(y1, x1, y2, x2):
                                rectangles.append((y1, x1, y2, x2))
    
    # Remove duplicates and return count
    unique_rectangles = set(rectangles)
    return len(unique_rectangles)

# Test grid
grid = """
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                  #########################################     
                                  #                                       #     
                                  #                                       #     
                                  #                                       #     
                                  #                                       #     
              ####################█##################################     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   #                                 #     #     
              #                   ##################################█######     
              #                                                     #           
              #                                                     #           
              #######################################################           
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
"""

print(find_rectangles(grid))