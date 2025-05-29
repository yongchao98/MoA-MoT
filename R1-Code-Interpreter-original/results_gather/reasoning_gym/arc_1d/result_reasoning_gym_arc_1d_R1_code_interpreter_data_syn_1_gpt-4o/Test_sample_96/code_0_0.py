def transform_grid(input_grid):
    # Move the first three numbers to the end
    return input_grid[3:] + input_grid[:3]

# Test input
test_input = [3, 7, 4, 4, 7, 9, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 2, 8, 1]
# Apply the transformation
output_grid = transform_grid(test_input)
print(output_grid)