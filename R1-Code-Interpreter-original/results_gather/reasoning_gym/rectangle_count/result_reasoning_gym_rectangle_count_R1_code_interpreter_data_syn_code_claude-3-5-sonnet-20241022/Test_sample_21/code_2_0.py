def count_rectangles(grid):
    # Convert grid to list of lists for faster access
    grid = [list(row) for row in grid]
    H, W = len(grid), len(grid[0])
    
    # Find all horizontal and vertical lines
    h_lines = []  # (y, x1, x2)
    v_lines = []  # (x, y1, y2)
    
    # Find horizontal lines
    for y in range(H):
        start = None
        for x in range(W):
            if grid[y][x] in '#█':
                if start is None:
                    start = x
            elif start is not None:
                h_lines.append((y, start, x-1))
                start = None
        if start is not None:
            h_lines.append((y, start, W-1))
    
    # Find vertical lines
    for x in range(W):
        start = None
        for y in range(H):
            if grid[y][x] in '#█':
                if start is None:
                    start = y
            elif start is not None:
                v_lines.append((x, start, y-1))
                start = None
        if start is not None:
            v_lines.append((x, start, H-1))
    
    # Count rectangles and overlaps
    rect_count = 0
    overlap_points = set()
    
    # For each horizontal line at top
    for y1, x1, x2 in h_lines:
        # Find matching bottom line
        for y2, x3, x4 in h_lines:
            if y2 <= y1:
                continue
            # Check if x-coordinates overlap
            left = max(x1, x3)
            right = min(x2, x4)
            if left > right:
                continue
            
            # Check if vertical lines exist
            left_exists = False
            right_exists = False
            for x, y3, y4 in v_lines:
                if y3 <= y1 and y4 >= y2:
                    if x == left:
                        left_exists = True
                    if x == right:
                        right_exists = True
            
            if left_exists and right_exists:
                rect_count += 1
                # Check for overlaps
                for y in range(y1, y2+1):
                    for x in range(left, right+1):
                        if grid[y][x] == '█':
                            overlap_points.add((y, x))
    
    print(rect_count + len(overlap_points))

# Read input
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

count_rectangles(grid)