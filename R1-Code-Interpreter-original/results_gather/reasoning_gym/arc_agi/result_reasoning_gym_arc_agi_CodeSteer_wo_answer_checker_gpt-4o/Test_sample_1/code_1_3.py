def transform_grid(input_grid):
    # Function to reflect a grid vertically
    def reflect_vertical(grid):
        return [row[::-1] for row in grid]
    
    # Function to reflect a grid horizontally
    def reflect_horizontal(grid):
        return grid[::-1]
    
    # Function to rotate a grid 180 degrees
    def rotate_180(grid):
        return [row[::-1] for row in grid[::-1]]
    
    # Initialize an empty 8x8 grid
    output_grid = [[0] * 8 for _ in range(8)]
    
    # Get the transformations
    top_left = input_grid
    top_right = reflect_vertical(input_grid)
    bottom_left = reflect_horizontal(input_grid)
    bottom_right = rotate_180(input_grid)
    
    # Fill the output grid with the quadrants
    for i in range(4):
        for j in range(4):
            output_grid[i][j] = top_left[i][j]
            output_grid[i][j+4] = top_right[i][j]
            output_grid[i+4][j] = bottom_left[i][j]
            output_grid[i+4][j+4] = bottom_right[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [4, 6, 4, 2],
    [2, 5, 2, 4],
    [6, 5, 6, 6],
    [6, 5, 6, 6]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))