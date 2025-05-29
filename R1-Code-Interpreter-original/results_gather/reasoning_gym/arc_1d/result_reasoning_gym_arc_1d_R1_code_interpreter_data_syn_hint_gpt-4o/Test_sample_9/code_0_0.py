def transform_grid(input_grid):
    output_grid = input_grid.copy()
    first_non_zero_found = False
    for i in range(len(input_grid)):
        if input_grid[i] != 0 and not first_non_zero_found:
            first_non_zero_found = True
        elif first_non_zero_found and input_grid[i] == 0:
            output_grid[i:] = [0] * (len(input_grid) - i)
            break
    return output_grid

test_input = [0, 5, 0, 5, 0, 0, 0, 0, 0, 0, 5, 0, 0, 5, 5, 5, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)