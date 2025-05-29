def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    visited = [[False] * width for _ in range(height)]
    count = 0
    
    def is_border(y, x):
        return 0 <= y < height and 0 <= x < width and grid[y][x] in ['#', '█']
    
    def find_rectangle(y, x):
        nonlocal count
        if not is_border(y, x) or visited[y][x]:
            return
        
        # Find right boundary
        right = x
        while right + 1 < width and is_border(y, right + 1):
            right += 1
            
        # Find bottom boundary
        bottom = y
        while bottom + 1 < height and is_border(bottom + 1, x):
            bottom += 1
            
        # Validate rectangle
        valid = True
        has_overlap = False
        
        # Check all borders
        for i in range(y, bottom + 1):
            if not is_border(i, x) or not is_border(i, right):
                valid = False
                break
            if grid[i][x] == '█' or grid[i][right] == '█':
                has_overlap = True
                
        if valid:
            for j in range(x, right + 1):
                if not is_border(y, j) or not is_border(bottom, j):
                    valid = False
                    break
                if grid[y][j] == '█' or grid[bottom][j] == '█':
                    has_overlap = True
        
        if valid:
            count += 1
            if has_overlap:
                count += 1
            
            # Mark as visited
            for i in range(y, bottom + 1):
                for j in range(x, right + 1):
                    visited[i][j] = True
    
    # Scan grid once
    for y in range(height):
        for x in range(width):
            if grid[y][x] in ['#', '█'] and not visited[y][x]:
                find_rectangle(y, x)
    
    return count

# Parse input
grid = []
while True:
    try:
        line = input()
        grid.append(list(line))
    except EOFError:
        break

# Get result
result = find_rectangles(grid)
print(f"<<<{result}>>>")