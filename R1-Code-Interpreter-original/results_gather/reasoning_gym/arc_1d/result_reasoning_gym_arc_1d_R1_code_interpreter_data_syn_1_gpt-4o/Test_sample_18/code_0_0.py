def transform_grid(input_grid):
    output_grid = input_grid[:]
    # Transform the first two '2's in the first row
    for i in range(2):
        if output_grid[i] == 2:
            output_grid[i] = 8
    # Transform the first two '2's in the second row
    for i in range(6, 8):
        if output_grid[i] == 2:
            output_grid[i] = 8
    return output_grid

input_grid = [2, 2, 2, 2, 2, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)