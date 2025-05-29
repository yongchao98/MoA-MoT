def transform_grid(input_grid):
    # Move the first element to the end of the list
    return input_grid[1:] + [input_grid[0]]

# Test input
test_input = [7, 9, 6, 4, 0, 0, 0, 8, 7, 1, 3, 5]
# Transform the test input
output_grid = transform_grid(test_input)
print(output_grid)