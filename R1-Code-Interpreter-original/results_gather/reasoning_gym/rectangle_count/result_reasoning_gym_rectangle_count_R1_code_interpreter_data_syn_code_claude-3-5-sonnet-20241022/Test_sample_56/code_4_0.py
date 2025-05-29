def count_rectangles(grid):
    count = 0
    height = len(grid)
    width = len(grid[0])
    
    # Pre-process to find all corner points
    corners = []
    for i in range(height):
        for j in range(width):
            if grid[i][j] in '#█':
                # Only store points that could be top-left corners
                if (i == 0 or grid[i-1][j] not in '#█') and (j == 0 or grid[i][j-1] not in '#█'):
                    corners.append((i, j))
    
    # For each potential top-left corner
    for top, left in corners:
        # Find immediate right edge
        right = left + 1
        while right < width and grid[top][right] in '#█':
            right += 1
        right -= 1
        
        if right == left:  # Not wide enough
            continue
            
        # Find immediate bottom edge
        bottom = top + 1
        while bottom < height and grid[bottom][left] in '#█':
            bottom += 1
        bottom -= 1
        
        if bottom == top:  # Not tall enough
            continue
            
        # Quick check for bottom-right corner
        if grid[bottom][right] not in '#█':
            continue
            
        # Verify rectangle is complete
        valid = True
        # Check right edge
        for i in range(top, bottom + 1):
            if grid[i][right] not in '#█':
                valid = False
                break
        
        if not valid:
            continue
            
        # Check bottom edge
        for j in range(left, right + 1):
            if grid[bottom][j] not in '#█':
                valid = False
                break
        
        if valid:
            count += 1
    
    return count

# Read grid efficiently
grid = []
try:
    while True:
        line = input()
        if line.strip():  # Skip empty lines
            grid.append(line)
except EOFError:
    pass

# Get count
result = count_rectangles(grid)
print(f"<<<{result}>>>")