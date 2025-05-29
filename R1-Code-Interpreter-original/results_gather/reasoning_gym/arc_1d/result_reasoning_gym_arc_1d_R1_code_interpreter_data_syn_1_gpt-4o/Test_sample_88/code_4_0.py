def transform_grid(input_grid):
    output_grid = input_grid[:]
    output_grid[1] = 0
    output_grid[4] = 0
    return output_grid

input_grid = [5, 0, 5, 5, 5, 0, 5, 5, 5, 5, 0, 5, 5, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)