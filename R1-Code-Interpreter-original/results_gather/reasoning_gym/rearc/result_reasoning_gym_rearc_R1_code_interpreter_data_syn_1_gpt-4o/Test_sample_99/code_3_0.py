def transform_grid(input_grid):
    # Find the position of the unique numbers (7 and 6 in this case)
    unique_positions = [(i, j) for i, row in enumerate(input_grid) for j, val in enumerate(row) if val in [7, 6]]
    
    # Determine the bounding box for cropping
    min_row = min(pos[0] for pos in unique_positions)
    max_row = max(pos[0] for pos in unique_positions)
    min_col = min(pos[1] for pos in unique_positions)
    max_col = max(pos[1] for pos in unique_positions)
    
    # Crop the grid
    cropped_grid = [row[min_col:max_col+1] for row in input_grid[min_row:max_row+1]]
    
    # Replace one of the unique numbers with 2
    replaced = False
    for i in range(len(cropped_grid)):
        for j in range(len(cropped_grid[i])):
            if cropped_grid[i][j] in [7, 6] and not replaced:
                cropped_grid[i][j] = 2
                replaced = True
    
    return cropped_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 7, 1, 1, 6]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))