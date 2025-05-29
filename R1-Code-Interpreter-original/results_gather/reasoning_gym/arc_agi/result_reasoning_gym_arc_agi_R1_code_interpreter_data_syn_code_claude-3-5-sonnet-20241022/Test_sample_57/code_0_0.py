def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create output grid with same values as input
    output = [row[:] for row in input_grid]
    
    # Add diagonal of 9's
    for i in range(rows-1):  # Stop before last row
        output[i][i] = 9
    
    # Fill last row with 1's except last column
    if rows > 1:  # Only if there's more than one row
        for j in range(cols-1):  # Don't change last column
            output[rows-1][j] = 1
    
    # Print the result in the required format
    result = []
    for row in output:
        result.append(' '.join(map(str, row)))
    print('\n'.join(result))

# Test input
test_input = [
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3]
]

transform_grid(test_input)