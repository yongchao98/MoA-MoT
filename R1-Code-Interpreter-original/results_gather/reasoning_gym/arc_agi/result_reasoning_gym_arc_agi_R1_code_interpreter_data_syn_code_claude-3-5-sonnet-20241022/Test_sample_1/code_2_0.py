def transform_grid(input_grid):
    # Create 8x8 output grid
    output = [[0]*8 for _ in range(8)]
    
    # Top-left quadrant: Original input
    for i in range(4):
        for j in range(4):
            output[i][j] = input_grid[i][j]
    
    # Top-right quadrant: Rotate clockwise
    for i in range(4):
        for j in range(4):
            output[i][j+4] = input_grid[3-j][i]
    
    # Bottom-left quadrant: Rotate counterclockwise
    for i in range(4):
        for j in range(4):
            output[i+4][j] = input_grid[j][3-i]
    
    # Bottom-right quadrant: Rotate 180 degrees
    for i in range(4):
        for j in range(4):
            output[i+4][j+4] = input_grid[3-i][3-j]
    
    return output

# Test input
input_grid = [
    [4, 6, 4, 2],
    [2, 5, 2, 4],
    [6, 5, 6, 6],
    [6, 5, 6, 6]
]

result = transform_grid(input_grid)

# Print the result in the required format
for row in result:
    print(' '.join(map(str, row)))