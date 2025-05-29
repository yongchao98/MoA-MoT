def transform_grid(input_grid):
    # Find the position of '2'
    pos_2 = input_grid.index(2)
    # Get the number immediately preceding '2'
    preceding_number = input_grid[pos_2 - 1]
    # Replace '2' and all subsequent numbers with the preceding number
    output_grid = input_grid[:pos_2] + [preceding_number] * (len(input_grid) - pos_2)
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0]
# Transform the test input
output = transform_grid(test_input)
print(output)