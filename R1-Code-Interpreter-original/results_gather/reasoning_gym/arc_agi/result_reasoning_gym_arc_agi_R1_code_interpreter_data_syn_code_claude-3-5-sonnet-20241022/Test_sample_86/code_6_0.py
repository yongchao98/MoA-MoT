def analyze_grid(input_grid):
    # Parse input grid
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Calculate output dimensions
    out_rows = rows // 3
    out_cols = cols // 2
    
    # Initialize output grid with 8s
    output = [[8] * out_cols for _ in range(out_rows)]
    
    def find_patterns_in_block(r, c):
        # Check only the exact 3x2 block starting at (r,c)
        patterns = []
        
        # Check horizontal patterns (must start within block)
        for row in range(r, r + 3):
            if row < rows and c + 2 < cols:
                if (grid[row][c] == 9 and 
                    grid[row][c + 1] == 9 and 
                    grid[row][c + 2] == 9):
                    patterns.append('H')
                    break  # Only count one horizontal pattern per block
        
        # Check vertical patterns (must be completely within block)
        for col in [c, c + 1]:
            if col < cols and r + 2 < rows:
                if (grid[r][col] == 9 and 
                    grid[r + 1][col] == 9 and 
                    grid[r + 2][col] == 9):
                    patterns.append('V')
                    break  # Only count one vertical pattern per block
        
        return patterns
    
    # Process each block
    for i in range(out_rows):
        for j in range(out_cols):
            r_start = i * 3
            c_start = j * 2
            
            patterns = find_patterns_in_block(r_start, c_start)
            
            # Apply patterns
            if 'H' in patterns:
                output[i][j] = 1
                if j + 1 < out_cols:
                    output[i][j + 1] = 3
            
            if 'V' in patterns:
                # Don't overwrite 1 from horizontal pattern
                if 'H' not in patterns:
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