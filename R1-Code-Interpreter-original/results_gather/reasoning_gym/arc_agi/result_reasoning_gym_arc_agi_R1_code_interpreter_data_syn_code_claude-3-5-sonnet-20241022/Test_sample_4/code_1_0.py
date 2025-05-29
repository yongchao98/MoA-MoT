def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    # Process each column
    for j in range(cols):
        # Get column values
        column = [input_grid[i][j] for i in range(rows)]
        
        # Check for multiple 1's or 5's
        for i in range(1, rows-1):  # Skip first and last rows
            if column[i] in [1, 5]:  # Only transform 1's and 5's
                if column.count(column[i]) > 1:  # If number appears multiple times
                    # Transform: 1 -> 5, 5 -> 1
                    output_grid[i][j] = 5 if column[i] == 1 else 1
    
    # Print the result
    for row in output_grid:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [1,1,8,5],
    [6,1,0,5],
    [6,0,2,4],
    [9,8,1,1]
]

transform_grid(test_input)