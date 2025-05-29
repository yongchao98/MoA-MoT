def find_rectangles(grid):
    rows, cols = len(grid), len(grid[0])
    
    # Find horizontal and vertical line segments
    h_lines = set()  # (start_y, start_x, end_x)
    v_lines = set()  # (start_x, start_y, end_y)
    
    # Find horizontal lines
    for i in range(rows):
        start = None
        for j in range(cols):
            if grid[i][j] in '#█':
                if start is None:
                    start = j
            elif start is not None:
                if j - start > 1:  # Minimum length of 2
                    h_lines.add((i, start, j-1))
                start = None
        if start is not None and cols - start > 1:
            h_lines.add((i, start, cols-1))
    
    # Find vertical lines
    for j in range(cols):
        start = None
        for i in range(rows):
            if grid[i][j] in '#█':
                if start is None:
                    start = i
            elif start is not None:
                if i - start > 1:  # Minimum length of 2
                    v_lines.add((j, start, i-1))
                start = None
        if start is not None and rows - start > 1:
            v_lines.add((j, start, rows-1))
    
    # Find rectangles by matching lines
    rectangles = set()
    for h1 in h_lines:  # top line
        y1, x1, x2 = h1
        for h2 in h_lines:  # bottom line
            y2, x3, x4 = h2
            if y2 <= y1 or x3 > x2 or x4 < x1:  # Skip invalid combinations
                continue
            
            # Check if vertical lines exist at both ends
            left_exists = False
            right_exists = False
            for v in v_lines:
                x, y_start, y_end = v
                if y_start <= y1 and y_end >= y2:
                    if x == x1:
                        left_exists = True
                    elif x == x2:
                        right_exists = True
                if left_exists and right_exists:
                    rectangles.add((y1, x1, y2, x2))
                    break
    
    return len(rectangles)

# Read input
grid = []
try:
    while True:
        line = input().rstrip('\n')
        if not line and grid:
            break
        grid.append(line)
except EOFError:
    pass

# Print result
print(f"<<<{find_rectangles(grid)}>>>")