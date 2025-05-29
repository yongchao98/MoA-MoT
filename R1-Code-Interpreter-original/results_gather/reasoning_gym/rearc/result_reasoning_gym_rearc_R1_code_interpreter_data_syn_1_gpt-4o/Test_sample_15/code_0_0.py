def transform_grid(input_grid):
    # Determine the size of the input grid
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    
    # Determine the size of the output grid
    if input_rows == 2 and input_cols == 2:
        output_size = 6
    elif input_rows == 3 and input_cols == 3:
        output_size = 21
    else:
        raise ValueError("Unexpected input grid size")
    
    # Create the output grid filled with zeros
    output_grid = [[0] * output_size for _ in range(output_size)]
    
    # Place the input grid in the top-left corner of the output grid
    for i in range(input_rows):
        for j in range(input_cols):
            output_grid[i][j] = input_grid[i][j]
    
    return output_grid

# Test input
input_grid = [
    [1, 0],
    [0, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))