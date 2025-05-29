def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is not None:
        # Move the first non-zero element to the start
        non_zero_elements = input_grid[first_non_zero_index:]
        zero_count = len(input_grid) - len(non_zero_elements)
        output_grid = non_zero_elements + [0] * zero_count
    else:
        # If all are zeros, return the same grid
        output_grid = input_grid
    
    return output_grid

# Test input
test_input = [8, 1, 2, 0, 0, 0, 0, 0, 0, 0, 8, 7, 2, 8, 4, 5, 1, 8, 2]
output = transform_grid(test_input)
print(output)