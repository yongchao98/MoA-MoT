def is_corner(grid, i, j):
    # Check if position is a corner ('#' or '█')
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    return grid[i][j] in ['#', '█']

def validate_rectangle(grid, top, left, bottom, right):
    # Check if the points form a valid rectangle
    if not all(is_corner(grid, top, j) for j in range(left, right + 1)):  # top edge
        return False
    if not all(is_corner(grid, bottom, j) for j in range(left, right + 1)):  # bottom edge
        return False
    if not all(is_corner(grid, i, left) for i in range(top, bottom + 1)):  # left edge
        return False
    if not all(is_corner(grid, i, right) for i in range(top, bottom + 1)):  # right edge
        return False
    return True

def count_rectangles(grid):
    height = len(grid)
    width = len(grid[0])
    rectangles = set()
    
    # Find all possible corners
    corners = [(i, j) for i in range(height) for j in range(width) 
              if grid[i][j] in ['#', '█']]
    
    # Check all possible rectangle combinations
    for top, left in corners:
        for bottom in range(top, height):
            for right in range(left, width):
                if is_corner(grid, bottom, right):
                    if validate_rectangle(grid, top, left, bottom, right):
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