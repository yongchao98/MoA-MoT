def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    def is_valid_char(c):
        return c in ['#', '█']
    
    def trace_rectangle(start_y, start_x):
        # Find width
        width_count = 0
        x = start_x
        while x < width and is_valid_char(grid[start_y][x]):
            width_count += 1
            x += 1
            
        # Find height
        height_count = 0
        y = start_y
        while y < height and is_valid_char(grid[y][start_x]):
            height_count += 1
            y += 1
            
        # Verify rectangle
        for y in range(start_y, start_y + height_count):
            for x in range(start_x, start_x + width_count):
                if not is_valid_char(grid[y][x]):
                    return None
                
        # Check if corners are present
        if (is_valid_char(grid[start_y][start_x]) and 
            is_valid_char(grid[start_y][start_x + width_count - 1]) and
            is_valid_char(grid[start_y + height_count - 1][start_x]) and
            is_valid_char(grid[start_y + height_count - 1][start_x + width_count - 1])):
            return (start_y, start_x, height_count, width_count)
        return None

    # Find all potential rectangles
    visited = set()
    for y in range(height):
        for x in range(width):
            if (y, x) not in visited and is_valid_char(grid[y][x]):
                rect = trace_rectangle(y, x)
                if rect:
                    rectangles.append(rect)
                    # Mark all points in this rectangle as visited
                    for ry in range(rect[0], rect[0] + rect[2]):
                        for rx in range(rect[1], rect[1] + rect[3]):
                            visited.add((ry, rx))

    # Count overlapping rectangles
    overlaps = set()
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                overlaps.add((y, x))
    
    # Each '█' character represents an overlap of 2 rectangles
    additional_rectangles = len(overlaps)
    
    return len(rectangles) + additional_rectangles

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

total_rectangles = find_rectangles(grid)
print(f"<<<{total_rectangles}>>>")