def is_corner(grid, x, y):
    if grid[y][x] not in '#█':
        return False
    
    # Count adjacent border characters
    adjacent_borders = 0
    if y > 0 and grid[y-1][x] in '#█':  # up
        adjacent_borders += 1
    if y < len(grid)-1 and grid[y+1][x] in '#█':  # down
        adjacent_borders += 1
    if x > 0 and grid[y][x-1] in '#█':  # left
        adjacent_borders += 1
    if x < len(grid[0])-1 and grid[y][x+1] in '#█':  # right
        adjacent_borders += 1
    
    return adjacent_borders == 2

def trace_rectangle(grid, x1, y1):
    height = len(grid)
    width = len(grid[0])
    
    # Find right edge
    x2 = x1
    while x2 < width and grid[y1][x2] in '#█':
        x2 += 1
    x2 -= 1
    
    # Find bottom edge
    y2 = y1
    while y2 < height and grid[y2][x1] in '#█':
        y2 += 1
    y2 -= 1
    
    # Verify rectangle
    if not (grid[y1][x2] in '#█' and grid[y2][x1] in '#█' and grid[y2][x2] in '#█'):
        return None
    
    # Check all borders are complete
    for x in range(x1, x2+1):
        if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
            return None
    for y in range(y1, y2+1):
        if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
            return None
            
    # Check interior (should be empty except for overlapping parts)
    for y in range(y1+1, y2):
        for x in range(x1+1, x2):
            if grid[y][x] == '#':
                return None
    
    return (x1, y1, x2, y2)

def find_rectangles(grid):
    # Convert grid to list of strings and remove empty lines
    lines = [line for line in grid.split('\n') if line]
    grid = [list(line) for line in lines]
    
    height = len(grid)
    width = len(grid[0])
    
    rectangles = set()
    
    # Find top-left corners and trace rectangles
    for y in range(height):
        for x in range(width):
            if grid[y][x] in '#█' and is_corner(grid, x, y):
                # Check if this could be a top-left corner
                if (y == 0 or grid[y-1][x] not in '#█') and (x == 0 or grid[y][x-1] not in '#█'):
                    rect = trace_rectangle(grid, x, y)
                    if rect:
                        rectangles.add(rect)
    
    return len(rectangles)

# The ASCII grid (same as before)
grid = """                                                                                
[... rest of the grid ...]
                                                                                """

print(find_rectangles(grid))