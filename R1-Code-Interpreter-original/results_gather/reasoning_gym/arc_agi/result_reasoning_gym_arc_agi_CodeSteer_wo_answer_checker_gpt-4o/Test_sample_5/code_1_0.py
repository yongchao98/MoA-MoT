def transform_grid(input_grid):
    # Determine the size of the input grid
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Identify rows and columns to extract based on non-7 numbers
    rows_to_extract = set()
    cols_to_extract = set()
    
    for i in range(num_rows):
        for j in range(num_cols):
            if input_grid[i][j] != 7:
                rows_to_extract.add(i)
                cols_to_extract.add(j)
    
    # Sort the rows and columns to maintain order
    rows_to_extract = sorted(rows_to_extract)
    cols_to_extract = sorted(cols_to_extract)
    
    # Form the output grid
    output_grid = []
    for i in rows_to_extract:
        output_row = [input_grid[i][j] for j in cols_to_extract]
        output_grid.append(output_row)
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 6, 7, 7, 7, 7, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 7, 1, 7, 7, 7, 7, 4, 2, 7, 4, 4, 7, 4, 9, 7, 4, 9, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 2, 7, 7, 7, 7, 7, 4, 2, 7, 4, 6, 7, 4, 6, 7, 4, 4, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 1, 7, 7, 7, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 4, 2, 7, 4, 6, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 4, 0, 7, 4, 4, 7, 4, 4, 7, 4, 0, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 9, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)

# Print the output grid in the required format
print("<<<output grid>>>")
for row in output_grid:
    print(" ".join(map(str, row)))