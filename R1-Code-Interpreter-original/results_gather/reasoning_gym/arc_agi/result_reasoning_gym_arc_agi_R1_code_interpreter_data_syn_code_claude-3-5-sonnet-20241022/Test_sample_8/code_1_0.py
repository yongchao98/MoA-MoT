def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def is_enclosed(i, j):
        if grid[i][j] == 8:
            return False
            
        # Check for vertical bounds
        left_bound = False
        right_bound = False
        left_pos = -1
        right_pos = -1
        
        # Find left bound
        for y in range(j-1, -1, -1):
            if grid[i][y] == 8:
                left_bound = True
                left_pos = y
                break
                
        # Find right bound
        for y in range(j+1, cols):
            if grid[i][y] == 8:
                right_bound = True
                right_pos = y
                break
        
        if not (left_bound and right_bound):
            return False
            
        # Check if there's a connecting horizontal 8 above or below
        # within the bounded region
        has_connection = False
        
        # Check above
        for x in range(i-1, -1, -1):
            if grid[x][j] == 8:
                break
            if grid[x][left_pos] == 8 and grid[x][right_pos] == 8:
                has_connection = True
                break
                
        if not has_connection:
            # Check below
            for x in range(i+1, rows):
                if grid[x][j] == 8:
                    break
                if grid[x][left_pos] == 8 and grid[x][right_pos] == 8:
                    has_connection = True
                    break
        
        return has_connection
    
    # Process each cell
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