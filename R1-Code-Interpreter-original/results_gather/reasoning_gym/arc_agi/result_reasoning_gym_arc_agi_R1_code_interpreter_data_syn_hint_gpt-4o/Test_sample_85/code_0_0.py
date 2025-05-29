def transform_grid(input_grid):
    # Find the positions of '1' and '9'
    pos_1 = None
    pos_9 = None
    for i, row in enumerate(input_grid):
        for j, value in enumerate(row):
            if value == 1:
                pos_1 = (i, j)
            elif value == 9:
                pos_9 = (i, j)
    
    # Create a copy of the input grid for the output
    output_grid = [row[:] for row in input_grid]
    
    # Draw a line of '8's from '1' to the column of '9'
    if pos_1 and pos_9:
        row_1, col_1 = pos_1
        _, col_9 = pos_9
        if col_1 < col_9:
            for col in range(col_1 + 1, col_9 + 1):
                output_grid[row_1][col] = 8
        elif col_1 > col_9:
            for col in range(col_1 - 1, col_9 - 1, -1):
                output_grid[row_1][col] = 8
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))