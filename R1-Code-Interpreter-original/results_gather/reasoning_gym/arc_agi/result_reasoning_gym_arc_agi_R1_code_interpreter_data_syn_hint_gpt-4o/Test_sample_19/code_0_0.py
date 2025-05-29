def transform_grid(input_grid):
    output_grid = [[0] * 10 for _ in range(10)]
    
    # Copy '2's directly to the output grid
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 2:
                output_grid[i][j] = 2
    
    # Transform '8's to new positions
    # Based on the examples, let's assume a transformation pattern
    # Example 1: (1,6) -> (4,5), (2,2) -> (4,7), (6,9) -> (6,8)
    # Example 2: (0,3) -> (2,5), (1,9) -> (2,8), (4,1) -> (5,6)
    # Let's assume a pattern: (i, j) -> (i+2, j+2) for '8's
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 8:
                new_i = (i + 2) % 10
                new_j = (j + 2) % 10
                output_grid[new_i][new_j] = 8
    
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