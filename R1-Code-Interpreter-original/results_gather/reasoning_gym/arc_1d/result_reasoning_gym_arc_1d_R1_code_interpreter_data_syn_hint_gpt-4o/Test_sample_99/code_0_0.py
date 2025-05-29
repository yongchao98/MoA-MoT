def transform_grid(input_grid):
    # Move the first element to the end and shift others to the left
    return input_grid[1:] + [input_grid[0]]

# Test input
input_grid = [0, 0, 0, 4, 0, 9, 0, 0, 3, 0, 0, 3, 0, 1, 0, 6, 0, 0, 0, 9, 0]
output_grid = transform_grid(input_grid)
print(output_grid)