def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Determine the central block to replace
    # For simplicity, assume the block is in the middle third of the grid
    start_row = rows // 3
    end_row = 2 * rows // 3
    start_col = cols // 3
    end_col = 2 * cols // 3
    
    # Create the output grid
    output_grid = []
    for i in range(rows):
        new_row = []
        for j in range(cols):
            if start_row <= i < end_row and start_col <= j < end_col:
                new_row.append(7)
            else:
                new_row.append(input_grid[i][j])
        output_grid.append(new_row)
    
    return output_grid

# Define the input grid
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

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))