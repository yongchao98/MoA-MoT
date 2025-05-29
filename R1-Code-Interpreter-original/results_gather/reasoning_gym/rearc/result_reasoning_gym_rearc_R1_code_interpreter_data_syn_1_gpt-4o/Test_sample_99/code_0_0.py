def transform_grid(input_grid):
    # Focus on the bottom part of the grid
    cropped_grid = input_grid[-1:]
    
    # Replace '7' with '2'
    output_grid = [[2 if cell == 7 else cell for cell in row] for row in cropped_grid]
    
    return output_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 7, 1, 1, 6]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))