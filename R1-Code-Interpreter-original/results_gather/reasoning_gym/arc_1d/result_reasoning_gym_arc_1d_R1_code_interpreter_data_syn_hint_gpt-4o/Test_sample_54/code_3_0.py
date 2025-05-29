def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If no non-zero element is found, return the input as output
    if first_non_zero_index is None:
        return input_grid
    
    # Find the sequence of non-zero elements
    non_zero_sequence = []
    for i in range(first_non_zero_index, len(input_grid)):
        if input_grid[i] == 0:
            break
        non_zero_sequence.append(input_grid[i])
    
    # Calculate the shift amount
    shift_amount = len(non_zero_sequence)
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    start_index = first_non_zero_index + shift_amount
    output_grid[start_index:start_index + len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)