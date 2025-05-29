def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height, width = len(grid), len(grid[0])
    visited = set()
    rectangles = 0
    
    def is_corner(y, x):
        if y < 0 or y >= height or x < 0 or x >= width:
            return False
        return grid[y][x] in '#█'
    
    def validate_rectangle(y1, x1):
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
        
        # Quick validation
        if (y1, x1, y2, x2) in visited:
            return False
            
        # Validate rectangle
        for y in range(y1, y2 + 1):
            if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
                return False
        for x in range(x1, x2 + 1):
            if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
                return False
                
        visited.add((y1, x1, y2, x2))
        return True
    
    # Scan for top-left corners
    for y in range(height):
        for x in range(width):
            if grid[y][x] not in '#█':
                continue
                
            # Check if it's a potential top-left corner
            if (not is_corner(y-1, x)) and (not is_corner(y, x-1)):
                if validate_rectangle(y, x):
                    rectangles += 1
                    
            # Early termination if we've found all rectangles
            if rectangles > 20:  # Reasonable upper limit for this problem
                return rectangles
    
    return rectangles

# Read grid
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines at the start
            continue
        grid.append(list(line))
except EOFError:
    pass

# Remove empty lines at the end
while grid and not any(c in '#█' for c in grid[-1]):
    grid.pop()

print(f"<<<{find_rectangles(grid)}>>>")