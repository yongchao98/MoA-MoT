def transform_grid(input_grid):
    output_grid = [0] * len(input_grid)
    first_non_zero_found = False
    for i in range(len(input_grid)):
        if input_grid[i] != 0 and not first_non_zero_found:
            first_non_zero_found = True
            output_grid[i] = input_grid[i]
        elif first_non_zero_found:
            if input_grid[i] != 0:
                output_grid[i] = input_grid[i]
            else:
                break
    return output_grid

test_input = [0, 5, 0, 5, 0, 0, 0, 0, 0, 0, 5, 0, 0, 5, 5, 5, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)