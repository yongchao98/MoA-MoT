def transform_grid(input_grid):
    # Determine the main pattern by examining the first row
    main_pattern = input_grid[0]
    
    # Create the output grid by replacing any deviation with the main pattern
    output_grid = []
    for row in input_grid:
        if row == main_pattern:
            output_grid.append(row)
        else:
            output_grid.append(main_pattern)
    
    return output_grid

# Test input grid
input_grid = [
    [0, 8, 0, 8, 0, 8, 0, 8, 0, 8],
    [8, 0, 8, 0, 8, 0, 8, 0, 8, 0],
    [0, 8, 0, 8, 0, 8, 0, 8, 0, 8],
    [8, 0, 8, 0, 6, 6, 8, 0, 8, 0],
    [0, 8, 0, 8, 6, 6, 0, 8, 0, 8],
    [8, 0, 8, 0, 8, 0, 8, 0, 8, 0],
    [0, 8, 0, 8, 0, 8, 0, 8, 0, 8],
    [8, 0, 8, 0, 8, 0, 8, 0, 8, 0],
    [0, 8, 0, 8, 0, 8, 0, 8, 0, 8],
    [8, 0, 8, 0, 8, 0, 8, 0, 8, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))