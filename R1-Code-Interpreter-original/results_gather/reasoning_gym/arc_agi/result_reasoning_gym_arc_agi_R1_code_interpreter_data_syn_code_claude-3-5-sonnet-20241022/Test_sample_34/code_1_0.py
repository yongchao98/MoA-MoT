def predict_output(input_grid):
    # Convert input to 2D list
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    
    # Create output grid (copy of input initially)
    output = [row[:] for row in grid]
    
    # First two rows remain unchanged
    # Third row becomes all 2's
    output[2] = [2, 2, 2, 2, 2]
    
    # For each column
    for col in range(5):
        # Check for 1,9 pattern vertically in input
        if grid[2][col] == 1 and grid[3][col] == 9:
            # Apply 2,9,1 pattern in output
            output[2][col] = 2
            output[3][col] = 9
            output[4][col] = 1
        else:
            # Keep 9's in row 4 where they are
            output[3][col] = grid[3][col]
            # Fill row 5 with 9's by default
            output[4][col] = 9
    
    # Print the predicted output
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = """2 2 2 2 2
2 2 2 2 2
1 2 1 2 2
9 2 9 2 2
9 9 9 9 9"""

predict_output(test_input)