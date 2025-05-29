def analyze_and_predict(input_grid):
    # Convert input to 2D list for easier manipulation
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    
    # Create output grid (copy of input initially)
    output = [row[:] for row in grid]
    
    # First two rows remain unchanged
    # Process rows 3-5
    
    # Reset row 3 to all 2's in output
    output[2] = [2] * 5
    
    # Keep 9's in row 4 where they are
    # For each column
    for col in range(5):
        # If there's a 1 in row 3 of input, put it in row 5 of output in same column
        if grid[2][col] == 1:
            output[4][col] = 1
            if grid[3][col] == 9:
                output[4][col] = 9
                output[4][col-1] = 9 if col > 0 else output[4][col]
                output[4][col+1] = 9 if col < 4 else output[4][col]
    
    # Print the predicted output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = """2 2 2 2 2
2 2 2 2 2
1 2 1 2 2
9 2 9 2 2
9 9 9 9 9"""

analyze_and_predict(test_input)