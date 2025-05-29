def transform_grid(input_grid):
    # Create output grid
    output = [[6 for _ in range(3)] for _ in range(3)]
    
    # Map the positions
    source_cols = [1, 3, 5]  # columns to extract from
    source_rows = [0, 2, 4]  # rows to extract from
    
    # Extract values from the specified positions
    for i, row in enumerate(source_rows):
        for j, col in enumerate(source_cols):
            output[i][j] = input_grid[row][col]
    
    return output

# Test input
test_input = [
    [6,4,6,4,6,0],
    [6,6,6,6,6,6],
    [6,6,6,6,6,0],
    [6,6,6,6,6,6],
    [6,1,6,1,6,7],
    [6,6,6,6,6,6]
]

result = transform_grid(test_input)
for row in result:
    print(row)