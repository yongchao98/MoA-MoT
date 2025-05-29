def transform_grid(input_grid):
    # Extract non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    # Calculate the number of leading zeros in the input
    leading_zeros = len(input_grid) - len(non_zero_elements)
    # Create the output grid by placing non-zero elements at the start
    output_grid = non_zero_elements + [0] * leading_zeros
    return output_grid

# Test input grid
test_input = [6, 3, 4, 2, 6, 9, 8, 1, 7, 0, 0, 0]
# Transform the test input grid
test_output = transform_grid(test_input)
print(test_output)