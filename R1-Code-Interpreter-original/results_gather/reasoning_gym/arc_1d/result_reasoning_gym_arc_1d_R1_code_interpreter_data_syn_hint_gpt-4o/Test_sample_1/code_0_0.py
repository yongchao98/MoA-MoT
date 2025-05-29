def transform_grid(input_grid):
    # Find the last two non-zero numbers
    last_two_non_zero = [num for num in input_grid if num != 0][-2:]
    # Create the output grid
    output_grid = last_two_non_zero + input_grid[:-2]
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 6, 4, 4, 7, 8, 8, 6, 9, 8, 5, 8, 2, 1, 8, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)