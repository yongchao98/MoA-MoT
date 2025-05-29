def transform_grid(input_grid):
    count_2 = input_grid.count(2)
    count_1 = input_grid.count(1)
    output_grid = [2] * count_2 + [1] * count_1 + [0] * (len(input_grid) - count_2 - count_1)
    return output_grid

# Test input grid
input_grid = [0, 1, 0, 0, 2, 1, 0, 0, 0, 1, 0, 1, 1, 2, 0]
output_grid = transform_grid(input_grid)
print(output_grid)