def is_valid_rectangle(grid, top, left, bottom, right):
    # Check top and bottom edges
    for x in range(left, right + 1):
        if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
            return False
    
    # Check left and right edges
    for y in range(top, bottom + 1):
        if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
            return False
    
    return True

def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Find all potential corners
    corners = []
    for y in range(height):
        for x in range(width):
            if grid[y][x] in ['#', '█']:
                corners.append((y, x))
    
    # Check all possible rectangle combinations
    for top, left in corners:
        for bottom in range(top, height):
            for right in range(left, width):
                if (grid[bottom][right] in ['#', '█'] and
                    grid[top][right] in ['#', '█'] and
                    grid[bottom][left] in ['#', '█']):
                    if is_valid_rectangle(grid, top, left, bottom, right):
                        rectangles.add((top, left, bottom, right))
    
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