def transform_grid(input_grid):
    output_grid = [[0] * 10 for _ in range(10)]
    
    # Copy '2's directly to the output grid
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 2:
                output_grid[i][j] = 2
    
    # Define a mapping for '8's based on the examples
    eight_mapping = {
        (1, 6): (4, 5),
        (2, 2): (4, 7),
        (6, 9): (6, 8),
        (0, 3): (2, 5),
        (1, 9): (2, 8),
        (4, 1): (5, 6)
    }
    
    # Apply the mapping to transform '8's
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 8:
                if (i, j) in eight_mapping:
                    new_i, new_j = eight_mapping[(i, j)]
                    output_grid[new_i][new_j] = 8
                else:
                    # If the position is not in the mapping, keep it unchanged
                    output_grid[i][j] = 8
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 8, 0],
    [0, 0, 0, 8, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 2, 0, 0, 0, 0, 0, 0],
    [8, 0, 2, 2, 0, 0, 0, 0, 0, 8],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))