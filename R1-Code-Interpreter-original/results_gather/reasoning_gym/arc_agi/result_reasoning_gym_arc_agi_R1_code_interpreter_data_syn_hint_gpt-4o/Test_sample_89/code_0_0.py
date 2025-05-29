def transform_grid(input_grid):
    # Determine the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Determine the middle section to transform
    # For simplicity, let's assume the middle section is 3 rows long, starting from row 10
    start_row = 10
    end_row = 13
    
    # Transform the middle section to rows of 7s
    for i in range(start_row, end_row):
        output_grid[i] = [7] * num_cols
    
    return output_grid

# Test input grid
input_grid = [
    [4, 4, 5, 5],
    [5, 4, 4, 4],
    [5, 4, 4, 5],
    [5, 4, 4, 4],
    [4, 5, 4, 4],
    [4, 5, 4, 5],
    [5, 4, 5, 5],
    [4, 4, 4, 5],
    [4, 5, 4, 5],
    [5, 5, 5, 4],
    [4, 5, 4, 5],
    [4, 4, 4, 5],
    [4, 5, 5, 5],
    [4, 4, 5, 5],
    [4, 4, 5, 5],
    [5, 5, 5, 4],
    [5, 4, 4, 5],
    [5, 4, 4, 4],
    [5, 4, 4, 5],
    [5, 5, 4, 4],
    [4, 5, 5, 5],
    [4, 5, 5, 5],
    [5, 5, 5, 5],
    [5, 5, 4, 5]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))