def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def should_fill(i, j):
        if grid[i][j] == 8:
            return False
            
        # Check horizontal enclosure
        left_8 = False
        left_pos = -1
        right_8 = False
        right_pos = -1
        
        # Find nearest 8's to the left and right
        for k in range(j-1, -1, -1):
            if grid[i][k] == 8:
                left_8 = True
                left_pos = k
                break
        for k in range(j+1, cols):
            if grid[i][k] == 8:
                right_8 = True
                right_pos = k
                break
                
        # Check vertical enclosure
        up_8 = False
        down_8 = False
        
        if left_8 and right_8:
            # Check for 8's above or below in the enclosed region
            for k in range(i-1, -1, -1):
                if grid[k][j] == 8 or (grid[k][left_pos] == 8 and grid[k][right_pos] == 8):
                    up_8 = True
                    break
            for k in range(i+1, rows):
                if grid[k][j] == 8 or (grid[k][left_pos] == 8 and grid[k][right_pos] == 8):
                    down_8 = True
                    break
            
            if up_8 or down_8:
                return True
        
        return False
    
    # Apply the filling rule
    for i in range(rows):
        for j in range(cols):
            if should_fill(i, j):
                output[i][j] = 0
    
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