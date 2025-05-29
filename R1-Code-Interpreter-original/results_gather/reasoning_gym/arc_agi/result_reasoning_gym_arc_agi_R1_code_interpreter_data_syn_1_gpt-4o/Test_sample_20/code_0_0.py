def expand_grid(input_grid):
    # Determine the size of the input grid
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    
    # Initialize the output grid with the appropriate size
    output_rows = input_rows * 3
    output_cols = input_cols * 3
    output_grid = [[2 for _ in range(output_cols)] for _ in range(output_rows)]
    
    # Fill the output grid based on the input grid
    for i in range(input_rows):
        for j in range(input_cols):
            # The center of the 3x3 block in the output grid
            output_grid[i*3 + 1][j*3 + 1] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [3, 2, 2, 3, 2, 5, 5, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [4, 2, 2, 3, 2, 5, 5, 5]
]

# Generate the output grid
output_grid = expand_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))