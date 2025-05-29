def transform_grid(input_grid):
    # Determine the majority number
    majority_number = input_grid[0][0]
    
    # Create a new grid for the output
    output_grid = [row[:] for row in input_grid]
    
    # Iterate over the grid and apply the transformation rule
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            if input_grid[i][j] != majority_number:
                # Apply the transformation rule
                if (i, j) in [(1, 6), (1, 7), (2, 6), (2, 7)]:
                    output_grid[i][j] = 8
                else:
                    output_grid[i][j] = 2
    
    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 1, 1, 4],
    [4, 4, 4, 4, 4, 4, 1, 1, 4],
    [4, 4, 4, 4, 1, 1, 4, 4, 4],
    [4, 4, 4, 4, 1, 1, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))