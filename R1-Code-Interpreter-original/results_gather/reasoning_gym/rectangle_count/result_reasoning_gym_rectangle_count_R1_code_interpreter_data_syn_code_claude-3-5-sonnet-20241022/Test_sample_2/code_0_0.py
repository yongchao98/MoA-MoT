def is_corner(grid, i, j):
    # Check if position is a corner (has perpendicular lines)
    if i < 0 or j < 0 or i >= len(grid) or j >= len(grid[0]):
        return False
    if grid[i][j] not in ['#', '█']:
        return False
    
    horizontal = False
    vertical = False
    
    # Check horizontal
    if j > 0 and grid[i][j-1] in ['#', '█']:
        horizontal = True
    if j < len(grid[0])-1 and grid[i][j+1] in ['#', '█']:
        horizontal = True
        
    # Check vertical
    if i > 0 and grid[i-1][j] in ['#', '█']:
        vertical = True
    if i < len(grid)-1 and grid[i+1][j] in ['#', '█']:
        vertical = True
        
    return horizontal and vertical

def find_rectangle(grid, start_i, start_j):
    # Find complete rectangle starting from top-left corner
    if not is_corner(grid, start_i, start_j):
        return None
    
    # Find width
    width = 1
    while start_j + width < len(grid[0]) and grid[start_i][start_j + width] in ['#', '█']:
        width += 1
    
    # Find height
    height = 1
    while start_i + height < len(grid) and grid[start_i + height][start_j] in ['#', '█']:
        height += 1
    
    # Verify rectangle
    for i in range(start_i, start_i + height):
        for j in range(start_j, start_j + width):
            if i in [start_i, start_i + height - 1] or j in [start_j, start_j + width - 1]:
                if grid[i][j] not in ['#', '█']:
                    return None
    
    return (start_i, start_j, height, width)

def count_rectangles(grid):
    rectangles = set()
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['#', '█']:
                rect = find_rectangle(grid, i, j)
                if rect:
                    rectangles.add(rect)
    return len(rectangles)

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(list(line))
    except EOFError:
        break

print(count_rectangles(grid))