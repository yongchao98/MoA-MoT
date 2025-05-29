def extract_subgrid(input_grid):
    # Find the bounding box of non-9 values
    min_row, max_row = len(input_grid), 0
    min_col, max_col = len(input_grid[0]), 0
    
    for r in range(len(input_grid)):
        for c in range(len(input_grid[0])):
            if input_grid[r][c] != 9:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)
    
    # Extract the 4x4 subgrid
    subgrid = []
    for r in range(min_row, min_row + 4):
        subgrid.append(input_grid[r][min_col:min_col + 4])
    
    return subgrid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 9, 9, 1, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 9, 9, 1, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 2, 2, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9]
]

# Extract and print the subgrid
output_grid = extract_subgrid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))