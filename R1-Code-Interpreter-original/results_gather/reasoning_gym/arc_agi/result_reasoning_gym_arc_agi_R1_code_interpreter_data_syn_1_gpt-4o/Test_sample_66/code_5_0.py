def extract_bottom_left(input_grid):
    # Determine the size of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the bottom-left corner where numbers are not 9
    start_row = rows
    start_col = cols
    
    # Start from the bottom-left and find the first non-9 element
    for r in range(rows-1, -1, -1):
        for c in range(cols):
            if input_grid[r][c] != 9:
                start_row = r
                start_col = c
                break
        if start_row != rows:
            break
    
    # Extract the subgrid
    output_grid = []
    for r in range(start_row, rows):
        row = []
        for c in range(start_col, cols):
            if input_grid[r][c] != 9:
                row.append(input_grid[r][c])
        if row:
            output_grid.append(row)
    
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 3, 0, 1, 3, 9],
    [9, 9, 9, 9, 9, 0, 3, 0, 0, 9],
    [9, 0, 0, 3, 9, 1, 0, 1, 3, 9],
    [9, 3, 0, 1, 9, 9, 9, 9, 9, 9],
    [9, 1, 0, 0, 9, 9, 9, 9, 9, 9],
    [9, 0, 3, 1, 9, 9, 0, 0, 0, 0],
    [9, 0, 0, 0, 9, 9, 3, 0, 1, 3],
    [9, 3, 0, 1, 9, 9, 0, 3, 0, 0],
    [9, 9, 9, 9, 9, 9, 0, 0, 0, 1]
]

output_grid = extract_bottom_left(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))