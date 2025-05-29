def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Transform each row
    for r in range(rows):
        non_six = [x for x in input_grid[r] if x != 6]
        input_grid[r] = [6] * (cols - len(non_six)) + non_six
    
    # Transform each column
    for c in range(cols):
        non_six = [input_grid[r][c] for r in range(rows) if input_grid[r][c] != 6]
        for r in range(rows):
            if r < rows - len(non_six):
                input_grid[r][c] = 6
            else:
                input_grid[r][c] = non_six[r - (rows - len(non_six))]
    
    return input_grid

# Test input grid
input_grid = [
    [6, 2, 6, 2, 2],
    [8, 6, 6, 8, 6],
    [6, 6, 9, 6, 6],
    [1, 6, 6, 1, 6],
    [3, 6, 6, 6, 6]
]

output_grid = transform_grid(input_grid)
print(output_grid)