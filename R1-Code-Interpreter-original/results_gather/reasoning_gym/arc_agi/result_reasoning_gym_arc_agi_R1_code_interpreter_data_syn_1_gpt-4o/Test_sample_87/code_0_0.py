def transform_grid(input_grid):
    # Find the unique column number
    unique_col_num = None
    for col in range(len(input_grid[0])):
        col_values = [input_grid[row][col] for row in range(len(input_grid))]
        if len(set(col_values)) == 1:
            unique_col_num = col_values[0]
            break

    # Find the unique row number
    unique_row_num = None
    for row in input_grid:
        if len(set(row)) == 1:
            unique_row_num = row[0]
            break

    # Create the output grid
    output_grid = [row[:] for row in input_grid]

    # Apply the transformation
    for row_index, row in enumerate(input_grid):
        if row[0] == unique_row_num:
            for col_index, value in enumerate(row):
                if value == unique_col_num:
                    output_grid[row_index][col_index] = unique_row_num

    return output_grid

# Test input grid
input_grid = [
    [2, 2, 7, 7, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 7, 7, 2, 2, 2, 2, 2, 2, 2],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    [2, 2, 7, 7, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 7, 7, 2, 2, 2, 2, 2, 2, 2]
]

output_grid = transform_grid(input_grid)
print(output_grid)