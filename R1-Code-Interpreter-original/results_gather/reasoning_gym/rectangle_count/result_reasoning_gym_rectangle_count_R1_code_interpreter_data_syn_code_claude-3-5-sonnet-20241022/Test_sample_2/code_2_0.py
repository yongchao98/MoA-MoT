def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    rows, cols = len(grid), len(grid[0])
    count = 0
    
    # Store coordinates of all '#' and '█' characters
    points = set()
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in ['#', '█']:
                points.add((i, j))
    
    # Function to check if a point is part of rectangle border
    def is_valid_rectangle(top, left, bottom, right):
        if (top, left) not in points or (top, right) not in points or \
           (bottom, left) not in points or (bottom, right) not in points:
            return False
            
        # Check top and bottom edges
        for j in range(left + 1, right):
            if (top, j) not in points or (bottom, j) not in points:
                return False
                
        # Check left and right edges
        for i in range(top + 1, bottom):
            if (i, left) not in points or (i, right) not in points:
                return False
        
        return True
    
    # Find top-left corners
    corners = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in ['#', '█']:
                # Check if it's a top-left corner
                if (i == 0 or grid[i-1][j] not in ['#', '█']) and \
                   (j == 0 or grid[i][j-1] not in ['#', '█']):
                    corners.append((i, j))
    
    # For each corner, find possible rectangles
    for top, left in corners:
        # Find right boundary
        right = left
        while right + 1 < cols and grid[top][right + 1] in ['#', '█']:
            right += 1
            
        # Find bottom boundary
        bottom = top
        while bottom + 1 < rows and grid[bottom + 1][left] in ['#', '█']:
            bottom += 1
            
        # Validate rectangle
        if is_valid_rectangle(top, left, bottom, right):
            count += 1
    
    return count

# Read the grid
grid = []
try:
    while True:
        line = input()
        if not line:
            break
        grid.append(list(line))
except EOFError:
    pass

print(find_rectangles(grid))