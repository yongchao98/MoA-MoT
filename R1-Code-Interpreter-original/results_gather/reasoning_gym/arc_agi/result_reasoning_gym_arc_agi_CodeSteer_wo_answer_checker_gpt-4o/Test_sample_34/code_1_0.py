def transform_grid(input_grid):
    # Initialize the output grid with the same dimensions as the input grid
    output_grid = [[2] * 5 for _ in range(5)]
    
    # Copy the fourth row from the input grid to the output grid
    output_grid[3] = input_grid[3]
    
    # Transform the fifth row
    for i in range(5):
        if input_grid[2][i] == 1:
            output_grid[4][i] = 1
        else:
            output_grid[4][i] = input_grid[4][i]
    
    return output_grid

# Test input grid
input_grid = [
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
    [1, 2, 1, 2, 2],
    [9, 2, 9, 2, 2],
    [9, 9, 9, 9, 9]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))