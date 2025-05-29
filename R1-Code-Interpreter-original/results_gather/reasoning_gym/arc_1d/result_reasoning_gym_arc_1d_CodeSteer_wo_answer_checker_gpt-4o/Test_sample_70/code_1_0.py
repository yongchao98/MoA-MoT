def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # Count the number of leading zeros before the first non-zero number
    leading_zeros = first_non_zero_index
    
    # Shift the sequence of non-zero numbers to the right by the number of leading zeros
    output_grid = [0] * len(input_grid)
    non_zero_sequence = [x for x in input_grid if x != 0]
    output_grid[leading_zeros:leading_zeros + len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
input_grid = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")