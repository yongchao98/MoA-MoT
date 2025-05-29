def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    def is_valid_char(c):
        return c in '#█'
    
    def scan_rectangle(top, left):
        # Find right boundary
        right = left + 1
        while right < width and is_valid_char(grid[top][right]):
            right += 1
        right -= 1
        
        if right == left:  # No horizontal line
            return False
        
        # Find bottom boundary
        bottom = top + 1
        while bottom < height and is_valid_char(grid[bottom][left]):
            bottom += 1
        bottom -= 1
        
        if bottom == top:  # No vertical line
            return False
            
        # Verify rectangle
        # Check right vertical line
        for y in range(top, bottom + 1):
            if not is_valid_char(grid[y][right]):
                return False
        
        # Check bottom horizontal line
        for x in range(left, right + 1):
            if not is_valid_char(grid[bottom][x]):
                return False
                
        return True
    
    # Scan for top-left corners
    for i in range(height - 1):
        for j in range(width - 1):
            if not is_valid_char(grid[i][j]):
                continue
                
            # Check if it's a potential top-left corner
            if (i == 0 or not is_valid_char(grid[i-1][j])) and \
               (j == 0 or not is_valid_char(grid[i][j-1])):
                if scan_rectangle(i, j):
                    count += 1
    
    return count

# Read the grid
grid = []
try:
    while True:
        line = input()
        if not line.strip():  # Skip empty lines at the end
            continue
        grid.append(list(line))
except EOFError:
    pass

print(find_rectangles(grid))