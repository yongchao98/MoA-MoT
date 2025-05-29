def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def check_u_shape(x, y):
        # Check for U-shape patterns (vertical)
        if y > 0 and y < cols - 1:
            # Check vertical U facing right
            if (x > 0 and x < rows - 1 and
                grid[x-1][y-1] == 8 and grid[x][y-1] == 8 and grid[x+1][y-1] == 8 and
                grid[x-1][y+1] == 8 and grid[x][y+1] == 8 and grid[x+1][y+1] == 8):
                return True
                
        # Check for horizontal U-shape
        if x > 0 and x < rows - 1:
            # Check horizontal U facing down
            if (y > 0 and y < cols - 1 and
                grid[x-1][y-1] == 8 and grid[x-1][y] == 8 and grid[x-1][y+1] == 8 and
                grid[x+1][y-1] == 8 and grid[x+1][y] == 8 and grid[x+1][y+1] == 8):
                return True
        return False

    def is_enclosed(x, y):
        if grid[x][y] == 8:
            return False
            
        # Check if point is between vertical 8's
        left_eight = False
        right_eight = False
        for j in range(y-1, -1, -1):
            if grid[x][j] == 8:
                left_eight = True
                break
        for j in range(y+1, cols):
            if grid[x][j] == 8:
                right_eight = True
                break
                
        # Check if point is between horizontal 8's
        up_eight = False
        down_eight = False
        for i in range(x-1, -1, -1):
            if grid[i][y] == 8:
                up_eight = True
                break
        for i in range(x+1, rows):
            if grid[i][y] == 8:
                down_eight = True
                break
        
        # Point should be filled if it's properly enclosed
        return ((left_eight and right_eight) or (up_eight and down_eight)) and check_u_shape(x, y)
    
    # First pass: mark enclosed areas
    for i in range(1, rows-1):
        for j in range(1, cols-1):
            if is_enclosed(i, j):
                output[i][j] = 0
    
    # Second pass: fill continuous regions
    changed = True
    while changed:
        changed = False
        for i in range(1, rows-1):
            for j in range(1, cols-1):
                if output[i][j] == 0:
                    for ni, nj in [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]:
                        if (output[ni][nj] != 8 and output[ni][nj] != 0 and 
                            grid[ni][nj] == grid[i][j]):
                            output[ni][nj] = 0
                            changed = True
    
    # Convert to string format
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
test_input = [
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,8,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,5,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,5,5,8,5,5,5,5],
    [5,5,5,5,5,5,5,8,8,8,8,5,5,5,5],
    [5,8,8,8,8,5,5,5,5,5,5,5,5,5,5],
    [5,8,5,5,8,5,5,5,8,8,8,8,5,5,5],
    [5,8,5,5,8,5,5,5,8,5,5,8,5,5,5],
    [5,8,5,8,8,5,5,5,8,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,8,8,8,8,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
]

find_enclosed_areas(test_input)