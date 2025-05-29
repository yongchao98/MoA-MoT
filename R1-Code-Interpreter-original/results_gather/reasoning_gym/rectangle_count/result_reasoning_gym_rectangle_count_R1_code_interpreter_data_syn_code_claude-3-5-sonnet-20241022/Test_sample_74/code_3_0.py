def count_rectangles(grid):
    # Convert to list and get dimensions
    grid = [list(line) for line in grid]
    rows, cols = len(grid), len(grid[0])
    rectangles = 0
    
    # Find rectangles by corners
    for r in range(rows):
        for c in range(cols):
            # Skip if not a potential corner
            if grid[r][c] not in '#█':
                continue
            
            # Quick check for minimum rectangle size
            if r >= rows-1 or c >= cols-1:
                continue
                
            # Find right edge
            right = c + 1
            while right < cols and grid[r][right] in '#█':
                right += 1
            if right - c < 2:  # Too narrow
                continue
                
            # Find bottom edge
            bottom = r + 1
            while bottom < rows and grid[bottom][c] in '#█':
                bottom += 1
            if bottom - r < 2:  # Too short
                continue
                
            # Verify rectangle
            valid = True
            has_overlap = False
            
            # Check corners first
            if (grid[r][right-1] not in '#█' or 
                grid[bottom-1][c] not in '#█' or 
                grid[bottom-1][right-1] not in '#█'):
                continue
            
            # Quick check for edges
            for i in range(r, bottom):
                if grid[i][c] not in '#█' or grid[i][right-1] not in '#█':
                    valid = False
                    break
                if grid[i][c] == '█' or grid[i][right-1] == '█':
                    has_overlap = True
            
            if not valid:
                continue
                
            for j in range(c, right):
                if grid[r][j] not in '#█' or grid[bottom-1][j] not in '#█':
                    valid = False
                    break
                if grid[r][j] == '█' or grid[bottom-1][j] == '█':
                    has_overlap = True
            
            if valid:
                rectangles += 1
                if has_overlap:
                    rectangles += 1
    
    print(rectangles)

# Read input
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines
            continue
        grid.append(line)
except EOFError:
    pass

count_rectangles(grid)