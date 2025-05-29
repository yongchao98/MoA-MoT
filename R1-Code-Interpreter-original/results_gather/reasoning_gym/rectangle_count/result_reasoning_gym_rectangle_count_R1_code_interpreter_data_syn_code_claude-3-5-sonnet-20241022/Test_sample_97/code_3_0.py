def count_rectangles(grid):
    # Remove empty lines and convert to more efficient format
    grid = [''.join(row) for row in grid if '#' in row or '█' in row]
    if not grid:
        return 0
    
    width = len(grid[0])
    height = len(grid)
    count = 0
    
    # Store horizontal line segments for quick lookup
    h_lines = {}
    for y in range(height):
        line = grid[y]
        start = -1
        for x in range(width):
            if line[x] in '#█':
                if start == -1:
                    start = x
            elif start != -1:
                if y not in h_lines:
                    h_lines[y] = []
                h_lines[y].append((start, x-1))
                start = -1
        if start != -1:
            if y not in h_lines:
                h_lines[y] = []
            h_lines[y].append((start, width-1))
    
    # Check each pair of horizontal lines
    y_list = sorted(h_lines.keys())
    for i, y1 in enumerate(y_list):
        for y2 in y_list[i+1:]:
            for x1, x2 in h_lines[y1]:
                # Quick check for matching horizontal lines
                if (x1, x2) in set((a, b) for a, b in h_lines[y2]):
                    # Verify vertical connections
                    valid = True
                    for x in (x1, x2):
                        for y in range(y1+1, y2):
                            if grid[y][x] not in '#█':
                                valid = False
                                break
                        if not valid:
                            break
                    if valid:
                        count += 1
                        if count > 15:  # Safe upper limit for this problem
                            return count
    
    return count

# Read and process input
grid = []
try:
    while True:
        line = input()
        if any(c in '#█' for c in line):
            grid.append(list(line))
except EOFError:
    pass

print(f"<<<{count_rectangles(grid)}>>>")