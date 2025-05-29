def transform_grid(input_grid):
    output_grid = []
    current_value = 0
    for value in input_grid:
        if value != 0:
            current_value = value
        output_grid.append(current_value)
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 9, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 7, 0, 0, 0, 5]
output_grid = transform_grid(input_grid)
print(output_grid)