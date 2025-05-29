def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    visited = [[False] * width for _ in range(height)]
    count = 0
    
    def is_border(y, x):
        return grid[y][x] in ['#', '█']
    
    def validate_rectangle(y1, x1):
        # Find the width of top edge
        x2 = x1
        while x2 < width and is_border(y1, x2):
            x2 += 1
        x2 -= 1
        
        if x2 == x1:  # Single point
            return False
        
        # Find the height of left edge
        y2 = y1
        while y2 < height and is_border(y2, x1):
            y2 += 1
        y2 -= 1
        
        if y2 == y1:  # Single point
            return False
            
        # Validate right edge
        for y in range(y1, y2 + 1):
            if not is_border(y, x2):
                return False
                
        # Validate bottom edge
        for x in range(x1, x2 + 1):
            if not is_border(y2, x):
                return False
        
        # Mark all border points as visited
        for y in range(y1, y2 + 1):
            visited[y][x1] = True
            visited[y][x2] = True
        for x in range(x1, x2 + 1):
            visited[y1][x] = True
            visited[y2][x] = True
            
        return True
    
    # Scan for unvisited border points that could be top-left corners
    for y in range(height - 1):
        for x in range(width - 1):
            if not visited[y][x] and is_border(y, x):
                # Check if this point is a top-left corner of a rectangle
                if validate_rectangle(y, x):
                    count += 1
    
    return count

# Read the grid efficiently
grid = []
try:
    while True:
        line = input()
        if line.strip():  # Skip empty lines
            grid.append(list(line))
except EOFError:
    pass

# Count rectangles
count = find_rectangles(grid)
print(f"<<<{count}>>>")