def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is None:
        return input_grid  # No non-zero element found, return the input as is
    
    # Find the sequence of non-zero numbers
    non_zero_sequence = []
    for i in range(first_non_zero_index, len(input_grid)):
        if input_grid[i] != 0:
            non_zero_sequence.append(input_grid[i])
        else:
            break
    
    # Shift the sequence to the left by one position
    output_grid = [0] * len(input_grid)
    output_grid[1:1+len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
test_input = [0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)