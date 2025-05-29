def count_rectangles_fast(grid):
    # Convert grid to list of strings for faster access
    grid = [list(row) for row in grid]
    height = len(grid)
    width = len(grid[0])
    
    # Count basic rectangles (looking only at corners)
    count = 0
    
    # Store coordinates of all '#' and '█' characters
    border_chars = set()
    overlap_chars = set()
    
    # Single pass to collect all border and overlap characters
    for y in range(height):
        for x in range(width):
            if grid[y][x] == '#':
                border_chars.add((y, x))
            elif grid[y][x] == '█':
                overlap_chars.add((y, x))
                border_chars.add((y, x))
    
    # Function to check if point exists in border_chars
    def is_border(y, x):
        return (y, x) in border_chars
    
    # Find rectangles by checking corners only
    corners = []
    for y, x in border_chars:
        # Check if it's a top-left corner
        if not is_border(y-1, x) and not is_border(y, x-1):
            corners.append((y, x))
    
    # For each corner, find matching rectangle
    for y, x in corners:
        # Find right edge
        x2 = x
        while x2 < width and is_border(y, x2):
            x2 += 1
        x2 -= 1
        
        # Find bottom edge
        y2 = y
        while y2 < height and is_border(y2, x):
            y2 += 1
        y2 -= 1
        
        # Quick validation of rectangle
        if x2 > x and y2 > y:
            if is_border(y2, x2) and is_border(y, x2):
                count += 1
    
    # Add overlapping rectangles
    count += len(overlap_chars)
    
    return count

# Read input efficiently
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

print(f"<<<{count_rectangles_fast(grid)}>>>")