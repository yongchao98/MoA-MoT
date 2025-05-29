def analyze_grid(input_grid):
    # Parse input grid
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Calculate output dimensions
    out_rows = rows // 3
    out_cols = cols // 2
    
    # Initialize output grid with 8s
    output = [[8] * out_cols for _ in range(out_rows)]
    
    # Process each output position
    for out_row in range(out_rows):
        for out_col in range(out_cols):
            # Get corresponding input block position
            in_row = out_row * 3
            in_col = out_col * 2
            
            # Check only first row of block for horizontal pattern
            if (in_col + 2 < cols and 
                grid[in_row][in_col] == 9 and 
                grid[in_row][in_col + 1] == 9 and 
                grid[in_row][in_col + 2] == 9):
                output[out_row][out_col] = 1
                if out_col + 1 < out_cols:
                    output[out_row][out_col + 1] = 3
            
            # Check middle row of block for horizontal pattern
            elif (in_col + 2 < cols and in_row + 1 < rows and
                  grid[in_row + 1][in_col] == 9 and 
                  grid[in_row + 1][in_col + 1] == 9 and 
                  grid[in_row + 1][in_col + 2] == 9):
                output[out_row][out_col] = 1
                if out_col + 1 < out_cols:
                    output[out_row][out_col + 1] = 3
            
            # Check bottom row of block for horizontal pattern
            elif (in_col + 2 < cols and in_row + 2 < rows and
                  grid[in_row + 2][in_col] == 9 and 
                  grid[in_row + 2][in_col + 1] == 9 and 
                  grid[in_row + 2][in_col + 2] == 9):
                output[out_row][out_col] = 1
                if out_col + 1 < out_cols:
                    output[out_row][out_col + 1] = 3
            
            # Check for vertical pattern in first column of block
            if (in_row + 2 < rows and
                grid[in_row][in_col] == 9 and 
                grid[in_row + 1][in_col] == 9 and 
                grid[in_row + 2][in_col] == 9):
                output[out_row][out_col] = 1
                if out_row + 1 < out_rows:
                    output[out_row + 1][out_col] = 3
            
            # Check for vertical pattern in second column of block
            if (in_row + 2 < rows and in_col + 1 < cols and
                grid[in_row][in_col + 1] == 9 and 
                grid[in_row + 1][in_col + 1] == 9 and 
                grid[in_row + 2][in_col + 1] == 9):
                output[out_row][out_col] = 1
                if out_row + 1 < out_rows:
                    output[out_row + 1][out_col] = 3
    
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