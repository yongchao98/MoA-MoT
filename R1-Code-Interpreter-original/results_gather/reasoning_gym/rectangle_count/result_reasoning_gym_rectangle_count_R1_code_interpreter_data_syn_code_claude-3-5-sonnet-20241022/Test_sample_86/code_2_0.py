def find_rectangles(grid_str):
    # Convert the grid string into a list of rows
    grid = [line for line in grid_str.strip().split('\n')]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def is_complete_rectangle(top, left, bottom, right):
        # Check if all corners exist
        if (grid[top][left] not in '#█' or grid[top][right] not in '#█' or
            grid[bottom][left] not in '#█' or grid[bottom][right] not in '#█'):
            return False
            
        # Check horizontal edges
        for x in range(left + 1, right):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
                
        # Check vertical edges
        for y in range(top + 1, bottom):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
                
        # Check if interior is empty or part of another rectangle
        for y in range(top + 1, bottom):
            for x in range(left + 1, right):
                if grid[y][x] not in ' #█':
                    return False
        
        return True

    # Find all rectangles
    rectangles = set()
    for y1 in range(height):
        for x1 in range(width):
            if grid[y1][x1] not in '#█':
                continue
            # Found potential top-left corner
            for y2 in range(y1 + 1, height):
                for x2 in range(x1 + 1, width):
                    if grid[y2][x2] not in '#█':
                        continue
                    # Found potential bottom-right corner
                    if is_complete_rectangle(y1, x1, y2, x2):
                        rectangles.add((y1, x1, y2, x2))

    return len(rectangles)

# Create the grid as a single string
grid = """                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
  #########################################                                     
  #                                       #                                     
  #                                       #                                     
  #########################################                                     
                                       #########                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                                
                                       #       #                             ###
                                       #########                             # #
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
                                                                                
                                                                                
                                                                                """

print(find_rectangles(grid))