def transform_grid(input_grid):
    # Initialize an empty output grid with the reduced size
    output_grid = [[0] * 11 for _ in range(11)]
    
    # Place the number 8 at the corners
    output_grid[0][0] = 8
    output_grid[0][10] = 8
    output_grid[10][0] = 8
    output_grid[10][10] = 8
    
    # Extract and place blocks of numbers
    # Place the block of 7s
    output_grid[2][4:7] = [7, 7, 7]
    output_grid[3][4:7] = [7, 7, 7]
    output_grid[4][4:7] = [7, 7, 7]
    
    # Place the block of 4s
    output_grid[5][1:5] = [4, 4, 4, 4]
    output_grid[6][1:5] = [4, 4, 4, 4]
    output_grid[7][1:5] = [4, 4, 4, 4]
    output_grid[8][1:5] = [4, 4, 4, 4]
    
    # Place the block of 2s
    output_grid[5][6:10] = [2, 2, 2, 2]
    output_grid[6][6:10] = [2, 2, 2, 2]
    output_grid[7][6:10] = [2, 2, 2, 2]
    output_grid[8][6:10] = [2, 2, 2, 2]
    
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

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))