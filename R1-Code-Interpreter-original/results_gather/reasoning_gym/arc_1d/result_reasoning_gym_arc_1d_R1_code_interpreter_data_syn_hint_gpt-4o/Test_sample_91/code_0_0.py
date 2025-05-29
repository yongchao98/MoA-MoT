def transform_grid(input_grid):
    # Extract non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    # Calculate the number of leading zeros
    leading_zeros = len(input_grid) - len(non_zero_elements)
    # Create the output grid
    output_grid = [0] * leading_zeros + non_zero_elements
    return output_grid

# Test input
test_input = [7, 2, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# Get the output
output = transform_grid(test_input)
print(output)