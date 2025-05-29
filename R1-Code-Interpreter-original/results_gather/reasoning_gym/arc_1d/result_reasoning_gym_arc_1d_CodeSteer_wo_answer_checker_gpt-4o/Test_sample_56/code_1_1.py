def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If there is no non-zero element, return the input as is
    if first_non_zero_index is None:
        return input_grid
    
    # Calculate the number of leading zeros in the input
    leading_zeros = first_non_zero_index
    
    # Create the output grid with the same length as the input
    output_grid = [0] * len(input_grid)
    
    # Shift the non-zero elements to the right by the number of leading zeros
    shift_amount = len(input_grid) - leading_zeros - (len(input_grid) - first_non_zero_index)
    for i in range(first_non_zero_index, len(input_grid)):
        if input_grid[i] != 0:
            output_grid[i + shift_amount] = input_grid[i]
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0]
test_output = transform_grid(test_input)

# Print the output in the required format
print("<<<" + " ".join(map(str, test_output)) + ">>>")