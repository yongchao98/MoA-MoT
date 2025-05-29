def find_enclosed_areas(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [row[:] for row in grid]
    
    def check_vertical_enclosure(i, j):
        # Check if there are 8's above and below
        up_8 = False
        down_8 = False
        
        # Look up
        for k in range(i-1, -1, -1):
            if grid[k][j] == 8:
                up_8 = True
                break
            if grid[k][j] != grid[i][j]:  # Different background number
                return False
                
        # Look down
        for k in range(i+1, rows):
            if grid[k][j] == 8:
                down_8 = True
                break
            if grid[k][j] != grid[i][j]:  # Different background number
                return False
                
        return up_8 and down_8
    
    def check_horizontal_enclosure(i, j):
        # Check if there are 8's to left and right
        left_8 = False
        right_8 = False
        
        # Look left
        for k in range(j-1, -1, -1):
            if grid[i][k] == 8:
                left_8 = True
                break
            if grid[i][k] != grid[i][j]:  # Different background number
                return False
                
        # Look right
        for k in range(j+1, cols):
            if grid[i][k] == 8:
                right_8 = True
                break
            if grid[i][k] != grid[i][j]:  # Different background number
                return False
                
        return left_8 and right_8
    
    # First pass: identify enclosed areas
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 8:  # Skip 8's
                if check_vertical_enclosure(i, j) or check_horizontal_enclosure(i, j):
                    output[i][j] = 0
    
    # Second pass: clean up isolated zeros
    for i in range(rows):
        for j in range(cols):
            if output[i][j] == 0:
                # Check if this zero is connected to other zeros or properly bounded
                has_connection = False
                for di, dj in [(-1,0), (1,0), (0,-1), (0,1)]:
                    ni, nj = i + di, j + dj
                    if (0 <= ni < rows and 0 <= nj < cols and 
                        (output[ni][nj] == 0 or output[ni][nj] == 8)):
                        has_connection = True
                        break
                if not has_connection:
                    output[i][j] = grid[i][j]
    
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