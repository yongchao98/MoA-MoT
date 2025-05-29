def is_rectangle(grid, top, left, bottom, right):
    # Check top and bottom edges
    for x in range(left, right + 1):
        if grid[top][x] not in '#█' or grid[bottom][x] not in '#█':
            return False
    
    # Check left and right edges
    for y in range(top, bottom + 1):
        if grid[y][left] not in '#█' or grid[y][right] not in '#█':
            return False
    
    return True

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = []
    
    # Find all potential rectangles
    for top in range(height):
        for left in range(width):
            if grid[top][left] not in '#█':
                continue
                
            # Find possible bottom-right corners
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if grid[bottom][right] not in '#█':
                        continue
                        
                    # Validate if it forms a rectangle
                    if is_rectangle(grid, top, left, bottom, right):
                        rectangles.append((top, left, bottom, right))
    
    return len(rectangles)

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(list(line))
    except EOFError:
        break

# Count rectangles
count = find_rectangles(grid)
print(f"<<<{count}>>>")