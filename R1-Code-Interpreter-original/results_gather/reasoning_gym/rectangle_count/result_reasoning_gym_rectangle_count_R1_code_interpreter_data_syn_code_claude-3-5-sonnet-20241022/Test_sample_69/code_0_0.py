def is_corner(grid, i, j):
    if grid[i][j] not in ['#', '█']:
        return False
    # Check if it's a corner by looking at adjacent cells
    horizontal = False
    vertical = False
    if j > 0 and grid[i][j-1] in ['#', '█']:
        horizontal = True
    if j < len(grid[i])-1 and grid[i][j+1] in ['#', '█']:
        horizontal = True
    if i > 0 and grid[i-1][j] in ['#', '█']:
        vertical = True
    if i < len(grid)-1 and grid[i+1][j] in ['#', '█']:
        vertical = True
    return horizontal and vertical

def is_valid_rectangle(grid, top, left, bottom, right):
    # Check if corners are marked
    if not all(grid[y][x] in ['#', '█'] for x, y in [(left, top), (right, top), (left, bottom), (right, bottom)]):
        return False
    
    # Check horizontal lines
    for x in range(left + 1, right):
        if grid[top][x] not in ['#', '█'] or grid[bottom][x] not in ['#', '█']:
            return False
    
    # Check vertical lines
    for y in range(top + 1, bottom):
        if grid[y][left] not in ['#', '█'] or grid[y][right] not in ['#', '█']:
            return False
    
    return True

def count_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    corners = []
    
    # Find all corners
    for i in range(height):
        for j in range(width):
            if is_corner(grid, i, j):
                corners.append((i, j))
    
    rectangles = set()
    # Check all possible corner pairs
    for i, (y1, x1) in enumerate(corners):
        for j, (y2, x2) in enumerate(corners[i+1:], i+1):
            # Only consider corners that could form a diagonal
            if x1 != x2 and y1 != y2:
                top = min(y1, y2)
                bottom = max(y1, y2)
                left = min(x1, x2)
                right = max(x1, x2)
                
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

# Count rectangles and print result
print(count_rectangles(grid))