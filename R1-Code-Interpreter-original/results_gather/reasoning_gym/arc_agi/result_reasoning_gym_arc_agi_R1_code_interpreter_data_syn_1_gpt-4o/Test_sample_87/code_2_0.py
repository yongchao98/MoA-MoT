def transform_grid(input_grid):
    # Find the unique column number
    unique_col_index = None
    for col in range(len(input_grid[0])):
        col_values = [input_grid[row][col] for row in range(len(input_grid))]
        if len(set(col_values)) == 1:
            unique_col_index = col
            break

    # Find the unique row number
    unique_row_index = None
    for row_index, row in enumerate(input_grid):
        if len(set(row)) == 1:
            unique_row_index = row_index
            break

    # Create the output grid
    output_grid = [row[:] for row in input_grid]

    # Apply the transformation
    if unique_col_index is not None and unique_row_index is not None:
        for row_index in range(len(output_grid)):
            if output_grid[row_index][unique_col_index] == input_grid[unique_row_index][unique_col_index]:
                output_grid[row_index][unique_col_index] = input_grid[unique_row_index][0]

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