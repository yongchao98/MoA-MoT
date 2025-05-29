def scan_grid(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Store horizontal and vertical line segments
    h_lines = []  # (start_x, end_x, y)
    v_lines = []  # (start_y, end_y, x)
    
    # Scan horizontal lines
    for y in range(height):
        start = None
        for x in range(width):
            if grid[y][x] in '#█':
                if start is None:
                    start = x
            elif start is not None:
                h_lines.append((start, x-1, y))
                start = None
        if start is not None:
            h_lines.append((start, width-1, y))
    
    # Scan vertical lines
    for x in range(width):
        start = None
        for y in range(height):
            if grid[y][x] in '#█':
                if start is None:
                    start = y
            elif start is not None:
                v_lines.append((start, y-1, x))
                start = None
        if start is not None:
            v_lines.append((start, height-1, x))
    
    # Find rectangles by matching horizontal and vertical lines
    for h1_idx, (h1_start, h1_end, h1_y) in enumerate(h_lines):
        for h2_idx, (h2_start, h2_end, h2_y) in enumerate(h_lines[h1_idx+1:], h1_idx+1):
            if h2_y <= h1_y:
                continue
            
            # Check if there are matching vertical lines
            left_found = False
            right_found = False
            
            for v_start, v_end, v_x in v_lines:
                if v_start <= h1_y and v_end >= h2_y:
                    if v_x == h1_start:
                        left_found = True
                    elif v_x == h1_end:
                        right_found = True
                    
                    if left_found and right_found:
                        count += 1
                        break
            
    return count

# Read the grid
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines
            continue
        grid.append(line)
except EOFError:
    pass

# Count rectangles
count = scan_grid(grid)
print(f"<<<{count}>>>")