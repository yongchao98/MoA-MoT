def transform_grid(input_grid):
    # Determine the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Find the pattern in the input grid
    pattern = [input_grid[i][-1] for i in range(num_rows)]
    
    # Create the output grid by repeating the pattern in the last column
    output_grid = []
    for i in range(num_rows):
        new_row = input_grid[i][:-1] + [pattern[i]]
        output_grid.append(new_row)
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 3, 3, 3],
    [3, 3, 3, 3, 2],
    [3, 3, 3, 3, 2],
    [3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))