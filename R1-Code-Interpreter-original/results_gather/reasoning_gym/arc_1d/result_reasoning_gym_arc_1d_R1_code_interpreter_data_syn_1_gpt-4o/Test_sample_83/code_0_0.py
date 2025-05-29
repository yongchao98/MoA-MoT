def transform_grid(input_grid):
    # Extract non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    # Sort in descending order
    sorted_elements = sorted(non_zero_elements, reverse=True)
    # Fill the output grid
    output_grid = sorted_elements + [0] * (len(input_grid) - len(sorted_elements))
    return output_grid

# Test input
input_grid = [2, 1, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 1, 1]
# Get the output
output_grid = transform_grid(input_grid)
print(output_grid)