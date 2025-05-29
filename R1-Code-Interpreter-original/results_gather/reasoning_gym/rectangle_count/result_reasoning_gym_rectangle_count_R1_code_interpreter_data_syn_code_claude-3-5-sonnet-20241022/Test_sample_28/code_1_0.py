def count_rectangles(grid):
    # Clean input and convert to list of lists for easier manipulation
    grid = [list(row) for row in grid if any(c in '#█' for c in row)]
    if not grid:
        return 0
    
    H, W = len(grid), len(grid[0])
    visited = set()
    
    def find_rectangle(y, x):
        if (y, x) in visited or grid[y][x] not in '#█':
            return None
        
        # Find right edge
        right = x
        while right < W and grid[y][right] in '#█':
            right += 1
        right -= 1
        
        # Find bottom edge
        bottom = y
        while bottom < H and grid[bottom][x] in '#█':
            bottom += 1
        bottom -= 1
        
        # Quick validation
        if not all(grid[y][right] in '#█' and grid[bottom][x] in '#█'):
            return None
        
        # Mark all border points as visited
        for i in range(y, bottom + 1):
            for j in range(x, right + 1):
                if grid[i][j] in '#█':
                    if i in (y, bottom) or j in (x, right):
                        visited.add((i, j))
        
        return True
    
    # Count base rectangles
    rectangles = 0
    for y in range(H):
        for x in range(W):
            if find_rectangle(y, x):
                rectangles += 1
    
    # Count overlaps (simplified)
    overlaps = sum(1 for y in range(H) for x in range(W) 
                  if grid[y][x] == '█' and 
                  y > 0 and y < H-1 and x > 0 and x < W-1 and
                  grid[y-1][x] in '#█' and grid[y+1][x] in '#█' and
                  grid[y][x-1] in '#█' and grid[y][x+1] in '#█')
    
    return rectangles + overlaps

# Read input
grid = []
try:
    while True:
        line = input()
        grid.append(line)
except EOFError:
    pass

print(count_rectangles(grid))