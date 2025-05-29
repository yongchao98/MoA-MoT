def shift_grid(input_grid, shift_by):
    output_grid = [0] * len(input_grid)
    for i, value in enumerate(input_grid):
        if value != 0:
            new_position = (i + shift_by) % len(input_grid)
            output_grid[new_position] = value
    return output_grid

input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0]
shift_by = 7
output_grid = shift_grid(input_grid, shift_by)
print(output_grid)