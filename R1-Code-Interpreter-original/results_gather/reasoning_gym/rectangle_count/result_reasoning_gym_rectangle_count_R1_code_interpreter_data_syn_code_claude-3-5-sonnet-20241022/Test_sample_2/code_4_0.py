def quick_count_rectangles(grid):
    rows = len(grid)
    cols = len(grid[0])
    found = set()
    
    # Only check first 100 rows and columns to ensure quick execution
    r_limit = min(100, rows)
    c_limit = min(100, cols)
    
    for i in range(r_limit):
        for j in range(c_limit):
            # Only process if we find a potential corner
            if grid[i][j] in ['#', '█']:
                # Quick check if it's likely a top-left corner
                if (i == 0 or grid[i-1][j] not in ['#', '█']) and (j == 0 or grid[i][j-1] not in ['#', '█']):
                    # Find immediate width and height
                    w = 1
                    while j + w < cols and grid[i][j + w] in ['#', '█']:
                        w += 1
                    
                    h = 1
                    while i + h < rows and grid[i + h][j] in ['#', '█']:
                        h += 1
                    
                    # Quick validation of opposite corner
                    if (i + h - 1 < rows and j + w - 1 < cols and 
                        grid[i + h - 1][j + w - 1] in ['#', '█']):
                        found.add((i, j, h, w))
    
    return len(found)

# Minimal grid reading
grid = []
try:
    while True:
        line = input()
        if not line.strip():
            continue
        grid.append(list(line))
except:
    pass

print(quick_count_rectangles(grid))