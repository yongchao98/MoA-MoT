def transform_grid(input_grid):
    # Find the first non-zero element
    start_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if start_index is None:
        return input_grid  # No non-zero elements to shift

    # Count the number of non-zero elements
    non_zero_count = sum(1 for x in input_grid if x != 0)

    # Create the output grid with zeros
    output_grid = [0] * len(input_grid)

    # Place the non-zero elements shifted one position to the left
    for i in range(non_zero_count):
        output_grid[start_index - 1 + i] = input_grid[start_index + i]

    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 9, 9, 9, 9]
# Transform the test input
test_output = transform_grid(test_input)
print(test_output)