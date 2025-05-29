def transform_grid(input_grid):
    # Find the first non-zero element
    start_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If no non-zero element is found, return the input grid as is
    if start_index is None:
        return input_grid
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[start_index:]
    
    # Create the output grid with zeros
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence shifted by 4 positions
    output_grid[4:4+len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
test_input = [7, 2, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)