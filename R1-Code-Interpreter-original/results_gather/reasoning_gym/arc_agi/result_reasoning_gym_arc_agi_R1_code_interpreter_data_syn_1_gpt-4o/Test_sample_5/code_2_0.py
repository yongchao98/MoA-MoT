def extract_output_grid(input_grid):
    # Initialize lists to store the indices of rows and columns to extract
    rows_to_extract = []
    cols_to_extract = []
    
    # Determine which rows to extract based on the presence of non-7 numbers
    for i, row in enumerate(input_grid):
        if any(cell != 7 for cell in row):
            rows_to_extract.append(i)
    
    # Determine which columns to extract based on the presence of non-7 numbers
    num_cols = len(input_grid[0])
    for j in range(num_cols):
        if any(input_grid[i][j] != 7 for i in range(len(input_grid))):
            cols_to_extract.append(j)
    
    # Create the output grid by extracting the identified rows and columns
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

output_grid = extract_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))