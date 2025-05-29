def extract_output_grid(input_grid):
    # Based on the examples, let's try to identify a pattern for row and column selection
    # Hypothetical pattern: select every second row and column starting from the first non-7 row/column
    rows_to_extract = [1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16]
    cols_to_extract = [1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16]
    
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