def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_border(char):
        return char in ['#', '█']
    
    def trace_rectangle(start_y, start_x):
        # Try to trace a rectangle starting from this point
        if not is_border(grid[start_y][start_x]):
            return None
        
        # Find right edge
        right_x = start_x
        while right_x < width and is_border(grid[start_y][right_x]):
            right_x += 1
        right_x -= 1
        
        # Find bottom edge
        bottom_y = start_y
        while bottom_y < height and is_border(grid[bottom_y][start_x]):
            bottom_y += 1
        bottom_y -= 1
        
        # Verify rectangle
        if not is_border(grid[start_y][right_x]) or not is_border(grid[bottom_y][start_x]):
            return None
            
        # Check if all edges exist
        for y in range(start_y, bottom_y + 1):
            if not (is_border(grid[y][start_x]) and is_border(grid[y][right_x])):
                return None
        for x in range(start_x, right_x + 1):
            if not (is_border(grid[start_y][x]) and is_border(grid[bottom_y][x])):
                return None
                
        return (start_y, start_x, bottom_y, right_x)
    
    # Find all rectangles
    rectangles = set()
    overlap_points = set()
    
    # First pass: find rectangles
    for y in range(height):
        for x in range(width):
            if is_border(grid[y][x]):
                # Only start tracing from points that could be top-left corners
                if (x == 0 or not is_border(grid[y][x-1])) and (y == 0 or not is_border(grid[y-1][x])):
                    rect = trace_rectangle(y, x)
                    if rect:
                        rectangles.add(rect)
    
    # Second pass: count overlaps
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                count = 0
                for (y1, x1, y2, x2) in rectangles:
                    if y1 <= y <= y2 and x1 <= x <= x2:
                        count += 1
                if count > 1:
                    overlap_points.add((y, x))
    
    return len(rectangles) + len(overlap_points)

# Process input
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines at the end
            continue
        grid.append(line)
except EOFError:
    pass

print(find_rectangles(grid))