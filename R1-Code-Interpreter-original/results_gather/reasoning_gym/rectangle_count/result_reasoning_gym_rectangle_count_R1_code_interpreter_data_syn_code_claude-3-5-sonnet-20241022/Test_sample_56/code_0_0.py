def is_corner(grid, i, j):
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    return grid[i][j] in '#█'

def find_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Convert grid to list of strings for easier processing
    grid = [list(row) for row in grid]
    
    # Find all possible rectangles
    for top in range(height):
        for left in range(width):
            if not is_corner(grid, top, left):
                continue
                
            # Find possible bottom right corners
            for bottom in range(top + 1, height):
                for right in range(left + 1, width):
                    if not (is_corner(grid, top, right) and 
                           is_corner(grid, bottom, left) and 
                           is_corner(grid, bottom, right)):
                        continue
                    
                    # Verify edges
                    valid = True
                    # Check vertical edges
                    for i in range(top + 1, bottom):
                        if not (grid[i][left] in '#█' and grid[i][right] in '#█'):
                            valid = False
                            break
                    
                    # Check horizontal edges
                    for j in range(left + 1, right):
                        if not (grid[top][j] in '#█' and grid[bottom][j] in '#█'):
                            valid = False
                            break
                    
                    if valid:
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