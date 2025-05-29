def is_corner(grid, i, j):
    # Check if position is a corner (has perpendicular lines)
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    
    char = grid[i][j]
    if char not in ['#', '█']:
        return False
        
    # Count directions where we have continuous lines
    directions = 0
    # Check right
    if j+1 < len(grid[0]) and grid[i][j+1] in ['#', '█']:
        directions += 1
    # Check left
    if j-1 >= 0 and grid[i][j-1] in ['#', '█']:
        directions += 1
    # Check down
    if i+1 < len(grid) and grid[i+1][j] in ['#', '█']:
        directions += 1
    # Check up
    if i-1 >= 0 and grid[i-1][j] in ['#', '█']:
        directions += 1
    
    return directions == 2 and (
        (j+1 < len(grid[0]) and grid[i][j+1] in ['#', '█'] and i+1 < len(grid) and grid[i+1][j] in ['#', '█']) or
        (j+1 < len(grid[0]) and grid[i][j+1] in ['#', '█'] and i-1 >= 0 and grid[i-1][j] in ['#', '█']) or
        (j-1 >= 0 and grid[i][j-1] in ['#', '█'] and i+1 < len(grid) and grid[i+1][j] in ['#', '█']) or
        (j-1 >= 0 and grid[i][j-1] in ['#', '█'] and i-1 >= 0 and grid[i-1][j] in ['#', '█'])
    )

def find_rectangle(grid, start_i, start_j):
    # Find the complete rectangle starting from top-left corner
    i, j = start_i, start_j
    width = height = 0
    
    # Find width
    while j < len(grid[0]) and grid[i][j] in ['#', '█']:
        width += 1
        j += 1
    
    # Find height
    j = start_j
    while i < len(grid) and grid[i][j] in ['#', '█']:
        height += 1
        i += 1
    
    # Verify rectangle
    for i in range(start_i, start_i + height):
        for j in range(start_j, start_j + width):
            if i == start_i or i == start_i + height - 1 or j == start_j or j == start_j + width - 1:
                if grid[i][j] not in ['#', '█']:
                    return False
    return True

def count_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    rectangles = 0
    rows, cols = len(grid), len(grid[0])
    
    # Find all corners that could be top-left corners of rectangles
    for i in range(rows):
        for j in range(cols):
            if is_corner(grid, i, j):
                # Check if this corner is top-left corner
                if (i == 0 or grid[i-1][j] not in ['#', '█']) and (j == 0 or grid[i][j-1] not in ['#', '█']):
                    if find_rectangle(grid, i, j):
                        rectangles += 1

    return rectangles

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(line)
    except EOFError:
        break

print(count_rectangles(grid))