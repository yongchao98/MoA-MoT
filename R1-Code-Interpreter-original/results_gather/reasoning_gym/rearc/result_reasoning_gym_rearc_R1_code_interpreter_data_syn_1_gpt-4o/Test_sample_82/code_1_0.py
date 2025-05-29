def extract_central_block(input_grid):
    # Determine the dimensions of the grid
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Find the central number by checking the middle rows
    central_number = None
    for row in input_grid:
        for number in row:
            if row.count(number) > 1 and number != row[0]:
                central_number = number
                break
        if central_number is not None:
            break
    
    # Find the start and end indices of the central block
    start_row, end_row = None, None
    start_col, end_col = None, None
    
    for i, row in enumerate(input_grid):
        if central_number in row:
            if start_row is None:
                start_row = i
            end_row = i
            if start_col is None:
                start_col = row.index(central_number)
            end_col = len(row) - row[::-1].index(central_number) - 1
    
    # Extract the central block
    output_grid = []
    for i in range(start_row, end_row + 1):
        output_grid.append(input_grid[i][start_col:end_col + 1])
    
    return output_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 2, 2, 2, 2, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
]

output_grid = extract_central_block(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))