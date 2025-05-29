def transform_grid(input_grid):
    # Find the first zero
    first_zero_index = input_grid.index(0)
    # Remove the first zero and shift elements
    output_grid = input_grid[:first_zero_index] + input_grid[first_zero_index+1:]
    # Append a zero at the end
    output_grid.append(0)
    return output_grid

# Test input
test_input = [3, 8, 4, 0, 0, 0, 7, 7, 1, 0, 5, 0, 0]
# Get the output
test_output = transform_grid(test_input)
print(test_output)