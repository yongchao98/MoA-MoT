def transform_input_to_output(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
    # Replace the first 'first_non_zero_index' elements with zeros
    output_grid = [0] * first_non_zero_index + input_grid[first_non_zero_index:]
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 8, 9, 4, 2, 5, 4, 7, 7, 2, 9, 7, 0, 0, 0]
# Generate output
test_output = transform_input_to_output(test_input)
print(test_output)