def transform_grid(input_grid):
    output_grid = input_grid.copy()
    first_zero_found = False
    for i in range(len(input_grid)):
        if input_grid[i] == 0:
            first_zero_found = True
        elif first_zero_found:
            output_grid[i] = 0
    return output_grid

# Test input
input_grid = [0, 5, 0, 5, 0, 0, 0, 0, 0, 0, 5, 0, 0, 5, 5, 5, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)