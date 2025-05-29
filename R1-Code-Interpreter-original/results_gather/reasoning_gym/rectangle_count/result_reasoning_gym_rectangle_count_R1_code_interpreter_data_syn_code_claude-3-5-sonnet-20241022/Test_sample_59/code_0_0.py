def find_rectangles(grid):
    if not grid:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Helper function to check if a rectangle exists at given position
    def is_rectangle(top, left, bottom, right):
        # Check top and bottom horizontal lines
        for x in range(left, right + 1):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        
        # Check left and right vertical lines
        for y in range(top, bottom + 1):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
        
        return True
    
    # Find all potential top-left corners
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█':
                # Find potential bottom-right corners
                for bottom in range(y + 1, height):
                    for right in range(x + 1, width):
                        if grid[bottom][right] in '#█':
                            if is_rectangle(y, x, bottom, right):
                                count += 1
                                
    return count

# Create the grid from the visual representation
grid = [
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                           ###  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           # #  ",
    "                                                                           ###  ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                          #################     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #               #     ",
    "                                                          #################     ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                ",
    "                                                                                "
]

print(find_rectangles(grid))