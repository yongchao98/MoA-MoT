def transform_grid(input_grid):
    output_grid = input_grid[:]
    change_to_eight = True
    for i in range(len(input_grid)):
        if input_grid[i] == 2 and change_to_eight:
            output_grid[i] = 8
        elif input_grid[i] == 0:
            change_to_eight = False
    return output_grid

# Test input grid
input_grid = [2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)