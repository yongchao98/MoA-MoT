def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If no non-zero number is found, return the input as is
    if first_non_zero_index is None:
        return input_grid
    
    # Find the sequence of repeated numbers
    sequence_start = first_non_zero_index
    sequence_end = sequence_start
    
    # Find the end of the sequence
    for i in range(sequence_start + 1, len(input_grid)):
        if input_grid[i] == input_grid[sequence_start]:
            sequence_end = i
        else:
            break
    
    # Calculate the length of the sequence
    sequence_length = sequence_end - sequence_start + 1
    
    # Create the output grid
    output_grid = input_grid[:]
    
    # Shift the sequence to the left, maintaining the position of any leading non-zero numbers
    shift_index = first_non_zero_index
    for i in range(sequence_length):
        output_grid[shift_index + i] = input_grid[sequence_start]
    
    # Replace the trailing numbers with 0s
    for i in range(shift_index + sequence_length, len(output_grid)):
        output_grid[i] = 0
    
    return output_grid

# Test input
test_input = [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0]
output = transform_grid(test_input)
print(output)