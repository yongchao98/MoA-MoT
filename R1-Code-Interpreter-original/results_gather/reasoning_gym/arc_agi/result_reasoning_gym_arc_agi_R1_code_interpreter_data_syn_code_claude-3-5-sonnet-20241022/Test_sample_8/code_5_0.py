def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def is_enclosed(i, j):
        if grid[i][j] == 8:
            return False
            
        # Check horizontal enclosure
        left_8 = False
        right_8 = False
        left_idx = j
        right_idx = j
        
        # Look left
        for k in range(j-1, -1, -1):
            if grid[i][k] == 8:
                left_8 = True
                left_idx = k
                break
                
        # Look right
        for k in range(j+1, cols):
            if grid[i][k] == 8:
                right_8 = True
                right_idx = k
                break
        
        # If horizontally enclosed, check if there are any 8's above or below
        # within the bounded region
        if left_8 and right_8:
            # Look up within bounded region
            for k in range(i-1, -1, -1):
                if any(grid[k][m] == 8 for m in range(left_idx, right_idx+1)):
                    return True
                    
            # Look down within bounded region
            for k in range(i+1, rows):
                if any(grid[k][m] == 8 for m in range(left_idx, right_idx+1)):
                    return True
        
        # Check vertical enclosure
        up_8 = False
        down_8 = False
        up_idx = i
        down_idx = i
        
        # Look up
        for k in range(i-1, -1, -1):
            if grid[k][j] == 8:
                up_8 = True
                up_idx = k
                break
                
        # Look down
        for k in range(i+1, rows):
            if grid[k][j] == 8:
                down_8 = True
                down_idx = k
                break
        
        # If vertically enclosed, check if there are any 8's to left or right
        # within the bounded region
        if up_8 and down_8:
            # Look left within bounded region
            for k in range(j-1, -1, -1):
                if any(grid[m][k] == 8 for m in range(up_idx, down_idx+1)):
                    return True
                    
            # Look right within bounded region
            for k in range(j+1, cols):
                if any(grid[m][k] == 8 for m in range(up_idx, down_idx+1)):
                    return True
        
        return False
    
    # First pass: mark enclosed areas
    for i in range(rows):
        for j in range(cols):
            if is_enclosed(i, j):
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