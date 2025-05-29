def transform_grid(input_grid):
    # Extract non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    # Count zeros
    zero_count = len(input_grid) - len(non_zero_elements)
    # Construct the output grid
    output_grid = non_zero_elements[:-1] + [0] * zero_count + [non_zero_elements[0]]
    return output_grid

# Test input
input_grid = [4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
output_grid = transform_grid(input_grid)
print(output_grid)