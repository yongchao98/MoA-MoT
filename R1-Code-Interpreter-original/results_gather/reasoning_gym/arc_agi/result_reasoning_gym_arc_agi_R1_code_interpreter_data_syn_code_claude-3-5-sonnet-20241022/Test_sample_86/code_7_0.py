def analyze_grid(input_grid):
    # Parse input grid
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Calculate output dimensions
    out_rows = rows // 3
    out_cols = cols // 2
    
    # Initialize output grid with 8s
    output = [[8] * out_cols for _ in range(out_rows)]
    
    # Process each block
    for i in range(out_rows):
        for j in range(out_cols):
            r = i * 3
            c = j * 2
            
            # Only check for horizontal '999' at the start of each block
            if c + 2 < cols and grid[r][c] == 9 and grid[r][c+1] == 9 and grid[r][c+2] == 9:
                output[i][j] = 1
                if j + 1 < out_cols:
                    output[i][j+1] = 3
            
            # Check for horizontal '999' in middle of block
            if c + 2 < cols and r + 1 < rows and grid[r+1][c] == 9 and grid[r+1][c+1] == 9 and grid[r+1][c+2] == 9:
                output[i][j] = 1
                if j + 1 < out_cols:
                    output[i][j+1] = 3
            
            # Check for horizontal '999' at bottom of block
            if c + 2 < cols and r + 2 < rows and grid[r+2][c] == 9 and grid[r+2][c+1] == 9 and grid[r+2][c+2] == 9:
                output[i][j] = 1
                if j + 1 < out_cols:
                    output[i][j+1] = 3
    
    # Print output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = """9 9 9 8 8 8 8 8 8
8 9 8 8 8 8 8 9 8
8 8 8 9 8 8 9 9 9
8 8 9 9 9 8 8 8 8
8 8 8 8 8 9 8 8 8
9 9 9 8 9 9 9 8 8
8 9 8 8 8 8 8 8 8
8 8 8 8 8 8 9 9 9
8 8 9 9 9 8 8 9 8
8 8 8 9 8 8 8 8 8"""

analyze_grid(test_input)