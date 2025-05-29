def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is None:
        return input_grid  # No non-zero numbers, return as is
    
    # Find the sequence of numbers after the first non-zero number
    sequence_start = first_non_zero_index + 1
    sequence_number = input_grid[sequence_start]
    
    # Find the end of the sequence
    sequence_end = sequence_start
    while sequence_end < len(input_grid) and input_grid[sequence_end] == sequence_number:
        sequence_end += 1
    
    # Create the output grid
    output_grid = input_grid[:]
    
    # Shift the sequence left by one position
    if sequence_end > sequence_start:
        output_grid[sequence_start:sequence_end-1] = input_grid[sequence_start+1:sequence_end]
        output_grid[sequence_end-1] = 0
    
    return output_grid

# Test input
test_input = [9, 9, 9, 9, 9, 9, 9, 9, 0, 0, 0, 2]
output = transform_grid(test_input)
print(output)