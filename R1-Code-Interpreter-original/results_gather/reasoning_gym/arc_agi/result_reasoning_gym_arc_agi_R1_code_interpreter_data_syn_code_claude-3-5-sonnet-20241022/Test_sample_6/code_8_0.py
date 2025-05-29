def transform_grid(grid):
    rows = len(grid)
    cols = len(grid[0])
    output = [[8] * cols for _ in range(rows)]
    
    # Process each column from bottom up
    for j in range(cols):
        # Check if column has pattern at bottom
        if grid[rows-1][j] in [6, 9]:
            # Keep all connected 6s and 9s in this column until gap
            for i in range(rows-1, -1, -1):
                if grid[i][j] in [6, 9]:
                    output[i][j] = grid[i][j]
                else:
                    # Stop at first gap (8)
                    break
    
    return output

# Test input
test_input = [
    [8,8,8,8,8,8,8,8,6,6,6,8],
    [8,6,6,6,6,6,6,8,6,6,6,8],
    [8,6,9,6,6,6,6,8,8,8,8,8],
    [8,6,6,6,6,9,6,8,6,6,6,6],
    [8,6,6,6,6,6,6,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,6,9,6],
    [8,8,8,8,6,6,6,8,6,6,6,6],
    [8,8,8,8,6,6,9,8,8,8,8,8],
    [8,8,8,8,6,6,6,8,8,8,8,8],
    [6,6,6,8,9,6,6,8,8,8,8,8],
    [9,6,6,8,8,8,8,8,6,9,6,9],
    [6,6,6,8,8,8,8,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,9,6,6]
]

# Print result
result = transform_grid(test_input)
for row in result:
    print(' '.join(map(str, row)))