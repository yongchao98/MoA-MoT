def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    
    # Convert grid to binary representation for easier processing
    visited = [[False] * width for _ in range(height)]
    
    def flood_fill(y, x):
        if y < 0 or y >= height or x < 0 or x >= width:
            return None
        if visited[y][x] or grid[y][x] not in '#█':
            return None
            
        # Find rectangle boundaries
        left, right, top, bottom = x, x, y, y
        
        # Find horizontal bounds
        i = x
        while i >= 0 and grid[y][i] in '#█':
            left = i
            i -= 1
        i = x
        while i < width and grid[y][i] in '#█':
            right = i
            i += 1
            
        # Find vertical bounds
        i = y
        while i >= 0 and grid[i][x] in '#█':
            top = i
            i -= 1
        i = y
        while i < height and grid[i][x] in '#█':
            bottom = i
            i += 1
            
        # Mark all points in the rectangle as visited
        for i in range(top, bottom + 1):
            for j in range(left, right + 1):
                if grid[i][j] in '#█':
                    visited[i][j] = True
        
        return (top, left, bottom, right)
    
    rectangles = set()
    overlap_count = 0
    
    # Find all rectangles
    for y in range(height):
        for x in range(width):
            if not visited[y][x] and grid[y][x] in '#█':
                rect = flood_fill(y, x)
                if rect:
                    rectangles.add(rect)
            if grid[y][x] == '█':
                overlap_count += 1
    
    # Calculate final count
    # Base rectangles + overlapping rectangles (overlap_count // 2)
    return len(rectangles) + overlap_count // 2

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

# Count rectangles
result = find_rectangles(grid)
print(result)