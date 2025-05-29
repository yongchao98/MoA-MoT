def transform_grid(input_grid):
    # Determine the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Initialize the number to extend
    number_to_extend = None
    
    # Check the last few columns of each row for a change
    for row in range(num_rows):
        for col in range(num_cols - 1, -1, -1):
            if input_grid[row][col] != input_grid[row][0]:
                number_to_extend = input_grid[row][col]
                break
        if number_to_extend is not None:
            break
    
    # Create the output grid by extending the number
    output_grid = []
    for row in range(num_rows):
        new_row = []
        for col in range(num_cols):
            if number_to_extend is not None:
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