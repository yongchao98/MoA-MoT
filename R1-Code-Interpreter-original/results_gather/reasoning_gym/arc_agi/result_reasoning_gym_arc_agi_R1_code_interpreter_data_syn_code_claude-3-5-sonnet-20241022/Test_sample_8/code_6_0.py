def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def find_compartment(i, j):
        # Find the boundaries of the current compartment
        left, right, up, down = j, j, i, i
        
        # Find horizontal boundaries
        for k in range(j-1, -1, -1):
            if grid[i][k] == 8:
                left = k
                break
            if k == 0 or grid[i][k] != grid[i][j]:
                left = None
                break
                
        for k in range(j+1, cols):
            if grid[i][k] == 8:
                right = k
                break
            if k == cols-1 or grid[i][k] != grid[i][j]:
                right = None
                break
        
        # Find vertical boundaries
        for k in range(i-1, -1, -1):
            if grid[k][j] == 8:
                up = k
                break
            if k == 0 or grid[k][j] != grid[i][j]:
                up = None
                break
                
        for k in range(i+1, rows):
            if grid[k][j] == 8:
                down = k
                break
            if k == rows-1 or grid[k][j] != grid[i][j]:
                down = None
                break
        
        return left, right, up, down
    
    def should_fill(i, j):
        if grid[i][j] == 8:
            return False
            
        left, right, up, down = find_compartment(i, j)
        
        # Check if we have both horizontal and vertical boundaries
        if (left is not None and right is not None and 
            any(grid[k][left] == 8 and grid[k][right] == 8 
                for k in [max(0, i-1), min(rows-1, i+1)])):
            return True
            
        if (up is not None and down is not None and 
            any(grid[up][k] == 8 and grid[down][k] == 8 
                for k in [max(0, j-1), min(cols-1, j+1)])):
            return True
            
        return False
    
    # Mark cells that should be filled
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