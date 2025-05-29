def transform_grid(input_grid):
    # Find the first non-zero sequence
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If no non-zero element is found, return the input as is
    if first_non_zero_index is None:
        return input_grid
    
    # Find the last non-zero element in the sequence
    last_non_zero_index = len(input_grid) - 1 - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero_index:last_non_zero_index + 1]
    
    # Extract the leading zeros
    leading_zeros = input_grid[:first_non_zero_index]
    
    # Extract the trailing zeros
    trailing_zeros = input_grid[last_non_zero_index + 1:]
    
    # Construct the output grid
    output_grid = leading_zeros + trailing_zeros + non_zero_sequence
    
    return output_grid

# Test input
test_input = [5, 5, 9, 7, 4, 1, 5, 6, 9, 2, 7, 8, 1, 7, 1, 6, 8, 9, 9, 8, 8, 5, 4, 7, 4, 9, 6, 7, 7]

# Apply the transformation
test_output = transform_grid(test_input)

# Print the output
print("<<<" + " ".join(map(str, test_output)) + ">>>")