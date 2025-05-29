def transform_grid(input_grid):
    output_grid = input_grid.copy()
    unique_numbers = set(input_grid) - {0}
    
    for number in unique_numbers:
        first_index = input_grid.index(number)
        last_index = len(input_grid) - 1 - input_grid[::-1].index(number)
        for i in range(first_index, last_index + 1):
            output_grid[i] = number
    
    return output_grid

input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 0, 9]
output_grid = transform_grid(input_grid)
print(output_grid)