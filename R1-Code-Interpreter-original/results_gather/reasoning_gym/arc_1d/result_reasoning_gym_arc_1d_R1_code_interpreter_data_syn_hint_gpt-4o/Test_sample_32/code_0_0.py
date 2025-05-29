def transform_grid(input_grid):
    output_grid = []
    first_sequence = True
    for num in input_grid:
        if num == 3:
            if first_sequence:
                output_grid.append(2)
            else:
                output_grid.append(1)
        else:
            output_grid.append(num)
            first_sequence = False
    return output_grid

input_grid = [3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)