def transform_grid(input_grid):
    # Initialize the output grid with 8s
    output_grid = [[8] * len(input_grid[0]) for _ in range(len(input_grid))]
    
    # Retain the pattern of 1s and 9s in the bottom-right section
    for i in range(7, len(input_grid)):
        for j in range(8, len(input_grid[0])):
            output_grid[i][j] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [8, 1, 9, 1, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 1, 1, 1, 8, 8, 8, 8, 1, 1, 1, 1, 8],
    [8, 1, 1, 1, 8, 8, 8, 8, 1, 1, 9, 1, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 1, 1, 1, 1, 8],
    [8, 8, 8, 9, 1, 1, 1, 8, 1, 1, 1, 1, 8],
    [8, 8, 8, 1, 1, 1, 1, 8, 1, 9, 1, 1, 8],
    [8, 8, 8, 1, 1, 9, 1, 8, 1, 1, 1, 1, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [1, 1, 1, 8, 8, 8, 1, 1, 1, 1, 8, 1, 1],
    [9, 1, 9, 8, 8, 8, 1, 1, 1, 1, 8, 1, 1],
    [1, 1, 1, 8, 8, 8, 1, 9, 1, 1, 8, 1, 1],
    [1, 1, 9, 8, 8, 8, 1, 1, 1, 1, 8, 8, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))