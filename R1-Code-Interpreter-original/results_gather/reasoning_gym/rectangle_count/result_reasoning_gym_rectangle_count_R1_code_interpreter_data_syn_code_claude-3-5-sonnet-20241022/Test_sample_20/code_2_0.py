def find_rectangles(grid):
    if not grid or not grid[0]:
        return 0
    
    height = len(grid)
    width = len(grid[0])
    count = 0
    
    def is_valid_char(ch):
        return ch in ['#', '█']
    
    def check_rectangle(top, left):
        # Find right boundary
        right = left + 1
        while right < width and is_valid_char(grid[top][right]):
            right += 1
        right -= 1
        
        if right == left:  # Not wide enough
            return False
        
        # Find bottom boundary
        bottom = top + 1
        while bottom < height and is_valid_char(grid[bottom][left]):
            bottom += 1
        bottom -= 1
        
        if bottom == top:  # Not tall enough
            return False
            
        # Validate rectangle
        for i in range(top, bottom + 1):
            if not (is_valid_char(grid[i][left]) and is_valid_char(grid[i][right])):
                return False
        
        for j in range(left, right + 1):
            if not (is_valid_char(grid[top][j]) and is_valid_char(grid[bottom][j])):
                return False
                
        return True
    
    # Only check positions that could be top-left corners
    for i in range(height - 1):
        for j in range(width - 1):
            if is_valid_char(grid[i][j]):
                # Check if it's a potential top-left corner
                if (i == 0 or not is_valid_char(grid[i-1][j])) and \
                   (j == 0 or not is_valid_char(grid[i][j-1])):
                    if check_rectangle(i, j):
                        count += 1
    
    return count

# Read the grid
grid = []
while True:
    try:
        line = input()
        grid.append(list(line))
    except EOFError:
        break

print(find_rectangles(grid))