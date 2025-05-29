def transform_grid(input_grid):
    # Determine the size of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create an empty output grid with the required size
    output_grid = []
    
    # Repeat the input grid to form the output grid
    for i in range(rows):
        for j in range(rows):
            new_row = []
            for k in range(cols):
                new_row.extend([input_grid[i][k]] * rows)
            output_grid.append(new_row)
    
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 3, 3],
    [9, 3, 9, 9],
    [9, 9, 9, 3],
    [3, 9, 3, 9]
]

# Transform the input grid to get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))