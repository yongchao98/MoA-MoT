def transform_grid(input_grid):
    # Create 8x8 output grid
    output = [[0]*8 for _ in range(8)]
    
    # First, fill in the left half (0-3 columns)
    # Top-left quadrant (0-3 rows, 0-3 columns): Original input
    for i in range(4):
        for j in range(4):
            output[i][j] = input_grid[i][j]
    
    # Bottom-left quadrant (4-7 rows, 0-3 columns): Vertical mirror of top-left
    for i in range(4):
        for j in range(4):
            output[7-i][j] = input_grid[i][j]
    
    # Then, fill in the right half (4-7 columns)
    # Top-right and bottom-right quadrants: Horizontal mirror of left half
    for i in range(8):
        for j in range(4):
            output[i][7-j] = output[i][j]
    
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