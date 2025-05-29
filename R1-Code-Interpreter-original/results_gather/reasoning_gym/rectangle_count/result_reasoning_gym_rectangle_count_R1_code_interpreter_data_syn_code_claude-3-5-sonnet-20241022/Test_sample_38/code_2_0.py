def find_rectangles(grid):
    # First, let's identify the main rectangles by their corners
    rectangles = []
    height = len(grid)
    width = len(grid[0])
    
    def is_corner(y, x):
        if not (0 <= y < height and 0 <= x < width):
            return False
        return grid[y][x] in '#█'
    
    def trace_rectangle(y1, x1):
        # Find the right edge
        x2 = x1
        while x2 < width and grid[y1][x2] in '#█':
            x2 += 1
        x2 -= 1
        
        # Find the bottom edge
        y2 = y1
        while y2 < height and grid[y2][x1] in '#█':
            y2 += 1
        y2 -= 1
        
        # Verify it's a complete rectangle
        if not all(grid[y1][x] in '#█' for x in range(x1, x2 + 1)): # top edge
            return None
        if not all(grid[y2][x] in '#█' for x in range(x1, x2 + 1)): # bottom edge
            return None
        if not all(grid[y][x1] in '#█' for y in range(y1, y2 + 1)): # left edge
            return None
        if not all(grid[y][x2] in '#█' for y in range(y1, y2 + 1)): # right edge
            return None
            
        return (y1, x1, y2, x2)
    
    # Find all rectangles
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Check if this could be a top-left corner
                if (y == 0 or grid[y-1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█'):
                    rect = trace_rectangle(y, x)
                    if rect:
                        rectangles.append(rect)
    
    # Remove duplicates
    unique_rectangles = set(rectangles)
    return len(unique_rectangles)

# Create the grid
grid = []
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                         #########                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                     ####█#######█###################           ")
grid.append("                                     #   #       #                  #           ")
grid.append("                                     #   #       #                  #           ")
grid.append("                                     #   #       #                  #           ")
grid.append("                                     #   #       #                  #           ")
grid.append("                                     #   #       #                  #           ")
grid.append("                                     ####█#######█###################           ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #       #                              ")
grid.append("                                         #########                              ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")
grid.append("                                                                                ")

print(find_rectangles(grid))