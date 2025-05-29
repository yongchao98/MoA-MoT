def find_rectangles(grid_str):
    # Convert the grid string into a list of rows
    grid = [line for line in grid_str.strip().split('\n')]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0
    
    def get_rectangle_dimensions(y, x):
        """Find the dimensions of a rectangle starting at (y, x)"""
        if grid[y][x] not in '#█':
            return None
        
        # Find width
        width = 0
        for dx in range(x, len(grid[0])):
            if grid[y][dx] not in '#█':
                break
            width += 1
            
        # Find height
        height = 0
        for dy in range(y, len(grid)):
            if grid[dy][x] not in '#█':
                break
            height += 1
            
        # Verify it's a complete rectangle
        for dy in range(y, y + height):
            if grid[dy][x + width - 1] not in '#█':
                return None
        for dx in range(x, x + width):
            if grid[y + height - 1][dx] not in '#█':
                return None
                
        # Check if interior is valid
        for dy in range(y + 1, y + height - 1):
            for dx in range(x + 1, x + width - 1):
                if grid[dy][dx] not in ' #█':
                    return None
        
        return (height, width)

    rectangles = set()
    # Scan for top-left corners
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Check if this is a top-left corner
                if (y == 0 or grid[y-1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█'):
                    dims = get_rectangle_dimensions(y, x)
                    if dims:
                        rectangles.add((y, x, dims[0], dims[1]))

    return len(rectangles)

# Create the grid as a single string (same as before)
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