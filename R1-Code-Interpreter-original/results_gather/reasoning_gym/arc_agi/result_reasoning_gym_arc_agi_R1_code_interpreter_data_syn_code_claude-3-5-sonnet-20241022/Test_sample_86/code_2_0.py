def analyze_grid(input_grid):
    # Parse input grid
    grid = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Calculate output dimensions
    out_rows = rows // 3
    out_cols = cols // 2
    
    # Initialize output grid with 8s
    output = [[8] * out_cols for _ in range(out_rows)]
    
    # Function to check for three consecutive 9s
    def check_three_nines(r, c, direction):
        if direction == 'horizontal':
            if c + 2 >= cols:
                return False
            return grid[r][c] == 9 and grid[r][c+1] == 9 and grid[r][c+2] == 9
        else:  # vertical
            if r + 2 >= rows:
                return False
            return grid[r][c] == 9 and grid[r+1][c] == 9 and grid[r+2][c] == 9
    
    # Process each column in output grid
    for out_col in range(out_cols):
        in_col = out_col * 2
        
        # Process each row in output grid
        for out_row in range(out_rows):
            in_row = out_row * 3
            
            # Check for horizontal pattern in first row of block
            if check_three_nines(in_row, in_col, 'horizontal'):
                output[out_row][out_col] = 1
                if out_col + 1 < out_cols:
                    output[out_row][out_col + 1] = 3
            
            # Check for vertical pattern
            if check_three_nines(in_row, in_col, 'vertical'):
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