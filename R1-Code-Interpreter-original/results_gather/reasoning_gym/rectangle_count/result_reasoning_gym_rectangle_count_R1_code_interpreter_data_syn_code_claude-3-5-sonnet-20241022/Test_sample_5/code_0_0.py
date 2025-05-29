def find_rectangles(grid):
    if not grid:
        return 0
    
    rows = len(grid)
    cols = len(grid[0])
    
    def is_valid_rect(top, left, bottom, right):
        # Check top and bottom edges
        for x in range(left, right + 1):
            if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
                return False
        # Check left and right edges
        for y in range(top, bottom + 1):
            if grid[y][left] not in '#█' or grid[y][right] not in '#█':
                return False
        return True
    
    rectangles = set()
    
    # Find all possible rectangles
    for top in range(rows):
        for left in range(cols):
            if grid[top][left] not in '#█':
                continue
                
            # Find possible bottom and right boundaries
            for bottom in range(top + 1, rows):
                if grid[bottom][left] not in '#█':
                    break
                for right in range(left + 1, cols):
                    if grid[top][right] not in '#█':
                        break
                    if is_valid_rect(top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
    return len(rectangles)

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

# Count rectangles
count = find_rectangles(grid)
print(f"<<<{count}>>>")