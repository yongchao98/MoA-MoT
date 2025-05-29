def transform_grid(input_grid):
    # Find the first and last non-zero elements
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    last_non_zero_index = len(input_grid) - 1 - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)
    
    # Calculate the number of trailing zeros in the input
    trailing_zeros = len(input_grid) - 1 - last_non_zero_index
    
    # Create the output grid with the same length as the input
    output_grid = [0] * len(input_grid)
    
    # Shift the non-zero elements to the right by the number of trailing zeros
    for i in range(first_non_zero_index, last_non_zero_index + 1):
        output_grid[i + trailing_zeros] = input_grid[i]
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0]
test_output = transform_grid(test_input)

# Print the output in the required format
print("<<<" + " ".join(map(str, test_output)) + ">>>")