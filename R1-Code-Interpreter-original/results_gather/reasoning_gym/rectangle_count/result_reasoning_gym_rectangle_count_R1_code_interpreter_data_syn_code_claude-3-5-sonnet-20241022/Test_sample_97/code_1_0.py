def count_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    count = 0
    seen = set()
    
    # Find top-left corners only
    for y in range(height-2):  # -2 to ensure minimum rectangle size
        for x in range(width-2):
            # Skip if not a potential corner
            if grid[y][x] not in '#█':
                continue
            
            # Quick check if it could be a top-left corner
            if y > 0 and grid[y-1][x] in '#█':
                continue
            if x > 0 and grid[y][x-1] in '#█':
                continue
            
            # Find right edge
            right = x + 1
            while right < width and grid[y][right] in '#█':
                right += 1
            right -= 1
            
            if right <= x:
                continue
                
            # Find bottom edge
            bottom = y + 1
            while bottom < height and grid[bottom][x] in '#█':
                bottom += 1
            bottom -= 1
            
            if bottom <= y:
                continue
            
            # Quick validation of bottom-right corner
            if grid[bottom][right] not in '#█':
                continue
                
            # Check if rectangle already found
            rect = (y, x, bottom, right)
            if rect in seen:
                continue
                
            # Minimal validation of edges
            valid = True
            for i in range(y, bottom+1):
                if grid[i][x] not in '#█' or grid[i][right] not in '#█':
                    valid = False
                    break
            if valid:
                for j in range(x, right+1):
                    if grid[y][j] not in '#█' or grid[bottom][j] not in '#█':
                        valid = False
                        break
            
            if valid:
                seen.add(rect)
                count += 1
                
            if count >= 15:  # Safe upper limit
                return count

    return count

# Minimal input processing
grid = []
try:
    while True:
        line = input().rstrip()
        if '#' in line or '█' in line:
            grid.append(line)
except EOFError:
    pass

print(f"<<<{count_rectangles(grid)}>>>")