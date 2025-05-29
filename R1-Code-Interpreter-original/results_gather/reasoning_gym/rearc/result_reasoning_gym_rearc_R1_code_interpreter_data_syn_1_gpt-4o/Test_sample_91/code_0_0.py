def transform_grid(input_grid):
    # Determine the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Find the number to extend
    # We assume the number to extend is in the last few rows or columns
    # Here, we check the last two rows for a change
    for col in range(num_cols):
        if input_grid[-1][col] != input_grid[-2][col]:
            number_to_extend = input_grid[-2][col]
            break
    
    # Create the output grid by extending the number
    output_grid = []
    for row in range(num_rows):
        new_row = []
        for col in range(num_cols):
            if row >= num_rows - 2 and col >= num_cols - 4:
                new_row.append(number_to_extend)
            else:
                new_row.append(input_grid[row][col])
        output_grid.append(new_row)
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 4, 4, 4, 4],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))