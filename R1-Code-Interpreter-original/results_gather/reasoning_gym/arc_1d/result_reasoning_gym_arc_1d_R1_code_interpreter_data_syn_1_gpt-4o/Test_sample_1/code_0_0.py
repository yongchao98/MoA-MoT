def transform_grid(input_grid):
    # Take the last two numbers
    last_two = input_grid[-2:]
    # Move them to the front and shift the rest
    output_grid = last_two + input_grid[:-2]
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 6, 4, 4, 7, 8, 8, 6, 9, 8, 5, 8, 2, 1, 8, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)