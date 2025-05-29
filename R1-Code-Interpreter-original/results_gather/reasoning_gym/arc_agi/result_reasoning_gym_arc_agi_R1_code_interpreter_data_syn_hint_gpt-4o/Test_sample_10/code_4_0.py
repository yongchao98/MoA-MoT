def transform_grid(input_grid):
    # Initialize the output grid with 8s
    output_grid = [[8] * len(input_grid[0]) for _ in range(len(input_grid))]
    
    # Determine the region to retain
    min_row, max_row = len(input_grid), 0
    min_col, max_col = len(input_grid[0]), 0
    
    for i in range(len(input_grid)):
        for j in range(len(input_grid[0])):
            if input_grid[i][j] in [1, 9]:
                min_row = min(min_row, i)
                max_row = max(max_row, i)
                min_col = min(min_col, j)
                max_col = max(max_col, j)
    
    # Retain the pattern of 1s and 9s in the determined region
    for i in range(min_row, max_row + 1):
        for j in range(min_col, max_col + 1):
            output_grid[i][j] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [8, 1, 9, 1, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 1, 1, 1, 8, 8, 8, 8, 1, 1, 1, 1, 8],
    [8, 1, 1, 1, 8, 8, 8, 8, 1, 1, 9, 1, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 1, 1, 1, 1, 8],
    [8, 8, 8, 9, 1, 1, 1, 8, 1, 1, 1, 1, 8],
    [8, 8, 8, 1, 1, 1, 1, 8, 1, 9, 1, 1, 8],
    [8, 8, 8, 1, 1, 9, 1, 8, 1, 1, 1, 1, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [1, 1, 1, 8, 8, 8, 1, 1, 1, 1, 8, 1, 1],
    [9, 1, 9, 8, 8, 8, 1, 1, 1, 1, 8, 1, 1],
    [1, 1, 1, 8, 8, 8, 1, 9, 1, 1, 8, 1, 1],
    [1, 1, 9, 8, 8, 8, 1, 1, 1, 1, 8, 8, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))