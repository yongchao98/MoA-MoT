def transform_grid(input_grid):
    output_grid = input_grid[:]
    first_zero_index = output_grid.index(0)
    for i in range(first_zero_index + 1, len(output_grid)):
        if output_grid[i] != 0:
            output_grid[first_zero_index] = output_grid[i]
            break
    output_grid[-1] = 0
    return output_grid

input_grid = [8, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 5, 2, 0, 0, 0, 3, 0, 8, 5, 0, 0, 0, 0, 0, 2, 8]
output_grid = transform_grid(input_grid)
print(output_grid)