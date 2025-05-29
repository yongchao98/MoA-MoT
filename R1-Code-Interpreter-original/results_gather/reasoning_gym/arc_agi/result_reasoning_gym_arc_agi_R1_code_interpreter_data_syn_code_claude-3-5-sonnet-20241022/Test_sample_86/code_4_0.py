def analyze_grid(input_grid):
    # Parse input grid
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Calculate output dimensions
    out_rows = rows // 3
    out_cols = cols // 2
    
    # Initialize output grid with 8s
    output = [[8] * out_cols for _ in range(out_rows)]
    
    def check_block_patterns(r_start, c_start):
        # Returns whether there's a horizontal '999' pattern in this block
        has_horizontal = False
        has_vertical = False
        
        # Check horizontal patterns (only in the exact block)
        for r in range(r_start, min(r_start + 3, rows)):
            for c in range(c_start, min(c_start + 2, cols - 2)):
                if (grid[r][c] == 9 and 
                    grid[r][c+1] == 9 and 
                    grid[r][c+2] == 9):
                    has_horizontal = True
        
        # Check vertical patterns (only in the exact block)
        for c in range(c_start, min(c_start + 2, cols)):
            for r in range(r_start, min(r_start + 1, rows - 2)):
                if (grid[r][c] == 9 and 
                    grid[r+1][c] == 9 and 
                    grid[r+2][c] == 9):
                    has_vertical = True
        
        return has_horizontal, has_vertical
    
    # Process each block
    for i in range(out_rows):
        for j in range(out_cols):
            r_start = i * 3
            c_start = j * 2
            
            h_pattern, v_pattern = check_block_patterns(r_start, c_start)
            
            # Apply patterns
            if h_pattern:
                output[i][j] = 1
                if j + 1 < out_cols:
                    output[i][j + 1] = 3
            
            if v_pattern:
                if output[i][j] != 1:  # Don't overwrite 1 from horizontal pattern
                    output[i][j] = 1
                if i + 1 < out_rows:
                    output[i + 1][j] = 3
    
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