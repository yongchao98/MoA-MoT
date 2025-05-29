def transform_grid(input_grid):
    # Initialize the output grid with zeros
    output_grid = [[0] * 10 for _ in range(11)]
    
    # Place '8' at the corners
    output_grid[0][0] = 8
    output_grid[0][9] = 8
    output_grid[10][0] = 8
    output_grid[10][9] = 8
    
    # Place '4' in the middle block
    for i in range(3, 7):
        for j in range(3, 7):
            output_grid[i][j] = 4
    
    # Place '2' in the lower part
    for i in range(7, 10):
        for j in range(3, 7):
            output_grid[i][j] = 2
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))