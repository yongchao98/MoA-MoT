def transform_grid(input_grid):
    n = len(input_grid)
    m = len(input_grid[0])
    
    # Count non-'7' numbers in each row and column
    row_counts = [sum(1 for x in row if x != 7) for row in input_grid]
    col_counts = [sum(1 for row in input_grid if row[j] != 7) for j in range(m)]
    
    # Determine whether to place '2' in the middle row or column
    if sum(row_counts) >= sum(col_counts):
        # Place '2' in the middle row
        middle_row = n // 2
        for j in range(m):
            input_grid[middle_row][j] = 2
    else:
        # Place '2' in the middle column
        middle_col = m // 2
        for i in range(n):
            input_grid[i][middle_col] = 2
    
    return input_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 4, 4, 4, 4, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 4, 4, 4, 4, 7, 7, 7, 7, 7],
    [7, 4, 4, 4, 4, 7, 7, 7, 7, 7],
    [7, 4, 4, 4, 4, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))