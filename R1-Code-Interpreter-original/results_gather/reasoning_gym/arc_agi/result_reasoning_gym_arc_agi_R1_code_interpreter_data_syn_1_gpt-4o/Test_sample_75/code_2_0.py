def find_subgrid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    min_row, max_row = rows, -1
    min_col, max_col = cols, -1
    
    # Find the bounds of the subgrid containing non-5 values
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] != 5:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)
    
    # Extract the subgrid
    subgrid = []
    for r in range(min_row, max_row + 1):
        subgrid.append(input_grid[r][min_col:max_col + 1])
    
    return subgrid

# Test input grid
input_grid = [
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 2, 2, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 0, 5, 5, 5, 5, 5, 1, 1, 5, 5],
    [5, 5, 0, 0, 5, 5, 5, 5, 5, 1, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
]

output_grid = find_subgrid(input_grid)
print(output_grid)