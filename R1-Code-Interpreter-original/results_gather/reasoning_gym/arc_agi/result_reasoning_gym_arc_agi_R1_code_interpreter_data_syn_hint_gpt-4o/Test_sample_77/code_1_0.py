def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Define the transformation rule based on the examples
    # In the examples, the central block of 2s is replaced by 9s
    # We will apply this transformation to the central part of the grid
    for i in range(5, 11):
        for j in range(5, 11):
            if input_grid[i][j] == 2:
                output_grid[i][j] = 9
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 1, 8, 0, 1, 5, 4, 4, 5, 1, 0, 8, 1, 3, 3],
    [3, 8, 8, 1, 1, 1, 4, 4, 4, 4, 1, 1, 1, 8, 8, 3],
    [1, 8, 4, 0, 5, 4, 4, 2, 2, 2, 4, 5, 0, 4, 8, 1],
    [8, 1, 0, 8, 4, 4, 8, 2, 2, 2, 4, 4, 8, 0, 1, 8],
    [0, 1, 2, 2, 2, 2, 4, 2, 2, 2, 3, 1, 4, 5, 1, 0],
    [1, 1, 2, 2, 2, 2, 4, 7, 7, 4, 3, 3, 4, 4, 1, 1],
    [5, 4, 2, 2, 2, 2, 7, 7, 7, 7, 4, 4, 8, 4, 4, 5],
    [4, 4, 2, 2, 2, 2, 7, 3, 3, 7, 7, 4, 8, 8, 4, 4],
    [4, 4, 8, 8, 4, 7, 7, 3, 3, 7, 7, 4, 8, 8, 4, 4],
    [5, 4, 4, 8, 4, 4, 7, 7, 7, 7, 4, 4, 8, 4, 4, 5],
    [1, 1, 4, 4, 3, 3, 4, 7, 7, 4, 3, 3, 4, 4, 1, 1],
    [0, 1, 5, 4, 1, 3, 4, 4, 4, 4, 3, 1, 4, 5, 1, 0],
    [8, 1, 0, 8, 4, 4, 8, 8, 8, 8, 4, 4, 8, 0, 1, 8],
    [1, 8, 4, 0, 5, 4, 4, 8, 8, 4, 4, 5, 0, 4, 8, 1],
    [3, 8, 8, 1, 1, 1, 4, 4, 4, 4, 1, 1, 1, 8, 8, 3],
    [3, 3, 1, 8, 0, 1, 5, 4, 4, 5, 1, 0, 8, 1, 3, 3]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))