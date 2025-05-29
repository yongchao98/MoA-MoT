def transform_grid(input_grid):
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    
    # Find the non-zero blocks
    for i in range(1, len(input_grid) - 1):
        for j in range(1, len(input_grid[i]) - 1):
            if input_grid[i][j] != 0:
                outer_value = input_grid[i][j]
                inner_value = input_grid[i+1][j+1]
                
                # Expand the outer block
                output_grid[i-1][j] = outer_value
                output_grid[i+1][j] = outer_value
                output_grid[i][j-1] = outer_value
                output_grid[i][j+1] = outer_value
                
                # Expand the inner block
                output_grid[i][j] = inner_value
                output_grid[i+1][j] = inner_value
                output_grid[i][j+1] = inner_value
                output_grid[i+1][j+1] = inner_value
                
                # Surround the inner block with outer block values
                output_grid[i-1][j-1] = outer_value
                output_grid[i-1][j+1] = outer_value
                output_grid[i+1][j-1] = outer_value
                output_grid[i+2][j] = outer_value
                output_grid[i][j+2] = outer_value
                output_grid[i+2][j+1] = outer_value
                output_grid[i+1][j+2] = outer_value
                output_grid[i+2][j+2] = outer_value
                
                return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 4, 4, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 4, 4, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))