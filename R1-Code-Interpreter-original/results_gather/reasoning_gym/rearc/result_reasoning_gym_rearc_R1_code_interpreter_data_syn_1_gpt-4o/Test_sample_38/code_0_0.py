def transform_grid(input_grid):
    # Determine the size of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create an empty output grid with the appropriate size
    output_grid = []
    
    # Repeat the input grid to form the output grid
    for i in range(2 * rows):
        output_row = []
        for j in range(2 * cols):
            # Determine the value to place in the output grid
            if (i % 2 == 1) and (j % 2 == 1):
                # Place 8 in the checkerboard pattern
                output_row.append(8)
            else:
                # Repeat the input grid values
                output_row.append(input_grid[i % rows][j % cols])
        output_grid.append(output_row)
    
    return output_grid

# Test input
input_grid = [
    [4, 5],
    [4, 4]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))