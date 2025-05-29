def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    # Find horizontal lines
    h_lines = []  # (y, start_x, end_x)
    for y in range(height):
        start = None
        for x in range(width):
            if grid[y][x] in ['#', '█']:
                if start is None:
                    start = x
            elif start is not None:
                if x - start > 1:  # Minimum width of 2
                    h_lines.append((y, start, x-1))
                start = None
        if start is not None and width - start > 1:
            h_lines.append((y, start, width-1))
    
    # Find rectangles using horizontal lines
    for i, (y1, x1_start, x1_end) in enumerate(h_lines):
        for y2, x2_start, x2_end in h_lines[i+1:]:
            if y2 <= y1:
                continue
            
            # Check if lines could form a rectangle
            if x1_start == x2_start and x1_end == x2_end:
                # Verify vertical lines
                valid = True
                has_overlap = False
                
                # Check left vertical line
                for y in range(y1, y2+1):
                    if grid[y][x1_start] not in ['#', '█']:
                        valid = False
                        break
                    if grid[y][x1_start] == '█':
                        has_overlap = True
                
                if not valid:
                    continue
                
                # Check right vertical line
                for y in range(y1, y2+1):
                    if grid[y][x1_end] not in ['#', '█']:
                        valid = False
                        break
                    if grid[y][x1_end] == '█':
                        has_overlap = True
                
                if valid:
                    # Check for overlaps in horizontal lines
                    for x in range(x1_start, x1_end+1):
                        if grid[y1][x] == '█' or grid[y2][x] == '█':
                            has_overlap = True
                            break
                    
                    count += 1
                    if has_overlap:
                        count += 1

    return count

# Parse the input grid
grid = []
while True:
    try:
        line = input()
        grid.append(list(line))
    except EOFError:
        break

# Count rectangles
result = find_rectangles(grid)
print(f"<<<{result}>>>")