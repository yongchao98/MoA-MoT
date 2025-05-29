def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    def is_valid_char(c):
        return c in '#█'
    
    def trace_rectangle(start_y, start_x):
        # Find width
        width_count = 0
        x = start_x
        while x < width and is_valid_char(grid[start_y][x]):
            width_count += 1
            x += 1
            
        if width_count < 2:
            return None
            
        # Find height
        height_count = 0
        y = start_y
        while y < height and is_valid_char(grid[y][start_x]):
            height_count += 1
            y += 1
            
        if height_count < 2:
            return None
            
        # Validate rectangle
        for y in range(start_y, start_y + height_count):
            for x in range(start_x, start_x + width_count):
                if not is_valid_char(grid[y][x]):
                    return None
                
        # Check if it's a valid rectangle (not just a line)
        if width_count < 2 or height_count < 2:
            return None
            
        return (start_y, start_x, height_count, width_count)
    
    rectangles = set()
    # Find all potential rectangles
    for y in range(height):
        for x in range(width):
            if is_valid_char(grid[y][x]):
                # Check if this could be top-left corner
                if (y == 0 or not is_valid_char(grid[y-1][x])) and \
                   (x == 0 or not is_valid_char(grid[y][x-1])):
                    rect = trace_rectangle(y, x)
                    if rect:
                        rectangles.add(rect)
    
    # Count overlapping rectangles
    overlaps = set()
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                count = 0
                for rect in rectangles:
                    ry, rx, rh, rw = rect
                    if ry <= y < ry + rh and rx <= x < rx + rw:
                        count += 1
                if count == 2:  # If point is shared by exactly 2 rectangles
                    overlaps.add((y, x))
    
    # Each overlap point represents one additional rectangle
    total_rectangles = len(rectangles) + len(overlaps)
    print(total_rectangles)

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

find_rectangles(grid)