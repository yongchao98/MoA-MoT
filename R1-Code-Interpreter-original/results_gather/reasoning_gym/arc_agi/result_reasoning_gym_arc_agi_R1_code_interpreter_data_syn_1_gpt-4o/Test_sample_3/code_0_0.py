def extract_output_grid(input_grid):
    # Define the rows and columns to extract
    rows_to_extract = [1, 2, 4]  # 0-indexed: 2nd, 3rd, and 5th rows
    cols_to_extract = [1, 3, 5]  # 0-indexed: 2nd, 4th, and 6th columns
    
    # Extract the output grid
    output_grid = []
    for row in rows_to_extract:
        output_row = [input_grid[row][col] for col in cols_to_extract]
        output_grid.append(output_row)
    
    return output_grid

# Test input grid
input_grid = [
    [6, 4, 6, 4, 6, 0],
    [6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 0],
    [6, 6, 6, 6, 6, 6],
    [6, 1, 6, 1, 6, 7],
    [6, 6, 6, 6, 6, 6]
]

# Get the output grid
output_grid = extract_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))