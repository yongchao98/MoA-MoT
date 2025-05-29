def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is None:
        return input_grid  # If no non-zero number is found, return the input as is
    
    # Extract the first non-zero number
    first_non_zero = input_grid[first_non_zero_index]
    
    # Find the sequence of non-zero numbers after the first non-zero number
    non_zero_sequence = [x for x in input_grid[first_non_zero_index + 1:] if x != 0]
    
    # Construct the output grid
    output_grid = [0] * len(input_grid)
    output_grid[0] = first_non_zero
    output_grid[1:1 + len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
input_grid = [1, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)