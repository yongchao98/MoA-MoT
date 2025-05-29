def find_corners(grid):
    height = len(grid)
    width = len(grid[0])
    corners = []
    
    # Find corners by checking for adjacent wall patterns
    for y in range(height):
        for x in range(width):
            if grid[y][x] not in '#█':
                continue
                
            # Check if it's a corner by looking at adjacent cells
            is_corner = False
            
            # Top-left corner
            if (x == 0 or grid[y][x-1] not in '#█') and (y == 0 or grid[y-1][x] not in '#█'):
                corners.append((x, y, 'TL'))
                continue
                
            # Top-right corner
            if (x == width-1 or grid[y][x+1] not in '#█') and (y == 0 or grid[y-1][x] not in '#█'):
                corners.append((x, y, 'TR'))
                continue
                
            # Bottom-left corner
            if (x == 0 or grid[y][x-1] not in '#█') and (y == height-1 or grid[y+1][x] not in '#█'):
                corners.append((x, y, 'BL'))
                continue
                
            # Bottom-right corner
            if (x == width-1 or grid[y][x+1] not in '#█') and (y == height-1 or grid[y+1][x] not in '#█'):
                corners.append((x, y, 'BR'))
    
    return corners

def verify_rectangle(grid, x1, y1, x2, y2):
    # Check top and bottom edges
    for x in range(x1, x2 + 1):
        if grid[y1][x] not in '#█' or grid[y2][x] not in '#█':
            return False
    
    # Check left and right edges
    for y in range(y1, y2 + 1):
        if grid[y][x1] not in '#█' or grid[y][x2] not in '#█':
            return False
    
    return True

def count_rectangles(grid):
    if not grid or not grid[0]:
        return 0
        
    corners = find_corners(grid)
    rectangles = set()
    
    # Group corners by type
    tl_corners = [(x, y) for x, y, t in corners if t == 'TL']
    br_corners = [(x, y) for x, y, t in corners if t == 'BR']
    
    # For each top-left corner, find matching bottom-right corners
    for tl_x, tl_y in tl_corners:
        for br_x, br_y in br_corners:
            if br_x > tl_x and br_y > tl_y:
                if verify_rectangle(grid, tl_x, tl_y, br_x, br_y):
                    rectangles.add((tl_x, tl_y, br_x, br_y))
    
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
print(count_rectangles(grid))