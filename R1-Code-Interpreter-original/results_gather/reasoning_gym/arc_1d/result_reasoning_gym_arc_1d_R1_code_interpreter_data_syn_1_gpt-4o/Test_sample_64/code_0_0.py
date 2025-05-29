def transform_grid(input_grid):
    # Find the index of '2'
    index_of_2 = input_grid.index(2)
    # Find the number that follows '2'
    replacement_number = input_grid[index_of_2 + 1]
    # Replace '2' with the replacement number
    input_grid[index_of_2] = replacement_number
    return input_grid

# Test input
test_input = [2, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# Transform the test input
output_grid = transform_grid(test_input)
print(output_grid)