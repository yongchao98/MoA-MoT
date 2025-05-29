def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Helper function to check if a point is part of rectangle border
    def is_border(ch):
        return ch in '#█'
    
    # Helper function to validate rectangle
    def validate_rectangle(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if not is_border(grid[top][x]) or not is_border(grid[bottom][x]):
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if not is_border(grid[y][left]) or not is_border(grid[y][right]):
                return False
        return True
    
    # Find all potential rectangles
    for y1 in range(height):
        for x1 in range(width):
            if is_border(grid[y1][x1]):
                # Found potential top-left corner
                for y2 in range(y1 + 1, height):
                    for x2 in range(x1 + 1, width):
                        if is_border(grid[y2][x2]):
                            # Check if this forms a valid rectangle
                            if validate_rectangle(y1, x1, y2, x2):
                                rectangles.append((y1, x1, y2, x2))
    
    # Process overlapping rectangles (marked by █)
    final_count = len(rectangles)
    
    # Add extra rectangle for each █ character (as it represents overlap)
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '█':
                final_count += 1
    
    # Since each overlapping point is counted twice (once in regular rectangles
    # and once in the █ counting), we need to divide the extra count by 2
    overlap_count = sum(1 for y in range(height) for x in range(width) if grid[y][x] == '█')
    final_count -= overlap_count // 2
    
    return final_count

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

# Count rectangles
result = find_rectangles(grid)
print(result)