def analyze_grid(input_grid):
    # Parse input grid
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Calculate output dimensions
    out_rows = rows // 3
    out_cols = cols // 2
    
    # Initialize output grid with 8s
    output = [[8] * out_cols for _ in range(out_rows)]
    
    # Process each 3x2 block
    for block_row in range(out_rows):
        for block_col in range(out_cols):
            # Get the current 3x2 block boundaries
            r_start = block_row * 3
            c_start = block_col * 2
            
            # Check for horizontal '999' patterns in all three rows of the block
            for r in range(r_start, r_start + 3):
                if c_start + 2 < cols:  # Ensure we don't go out of bounds
                    if grid[r][c_start] == 9 and grid[r][c_start + 1] == 9 and grid[r][c_start + 2] == 9:
                        output[block_row][block_col] = 1
                        if block_col + 1 < out_cols:
                            output[block_row][block_col + 1] = 3
            
            # Check for vertical '999' patterns in both columns of the block
            for c in range(c_start, min(c_start + 2, cols)):
                if r_start + 2 < rows:  # Ensure we don't go out of bounds
                    if grid[r_start][c] == 9 and grid[r_start + 1][c] == 9 and grid[r_start + 2][c] == 9:
                        output[block_row][block_col] = 1
                        if block_row + 1 < out_rows:
                            output[block_row + 1][block_col] = 3

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