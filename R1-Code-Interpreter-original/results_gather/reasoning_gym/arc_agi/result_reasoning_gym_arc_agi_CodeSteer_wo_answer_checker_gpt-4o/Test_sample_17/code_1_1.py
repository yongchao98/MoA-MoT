def transform_grid(input_grid):
    # Initialize an empty grid with the same dimensions
    output_grid = [[4 for _ in range(10)] for _ in range(10)]
    
    # Define the transformation rules
    transformation = {8: 0, 9: 9, 0: 5, 2: 9, 5: 2}
    
    # Iterate over the input grid
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] in transformation:
                # Shift the cluster diagonally downwards and to the right
                if i + 1 < 10 and j + 1 < 10:
                    output_grid[i + 1][j + 1] = transformation[input_grid[i][j]]
    
    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 8, 8, 4, 4, 4, 8, 8, 4],
    [4, 8, 8, 8, 8, 4, 8, 8, 8, 8],
    [4, 4, 8, 4, 4, 4, 4, 8, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 2, 2, 4],
    [4, 4, 8, 8, 4, 4, 1, 1, 2, 1],
    [4, 8, 8, 8, 8, 4, 4, 1, 4, 4],
    [4, 4, 8, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))