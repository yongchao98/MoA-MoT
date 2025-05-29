def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # Find the sequence of repeated numbers
    repeated_number = None
    repeated_start = None
    for i in range(first_non_zero_index + 1, len(input_grid)):
        if input_grid[i] != 0:
            repeated_number = input_grid[i]
            repeated_start = i
            break
    
    # Find the end of the repeated sequence
    repeated_end = repeated_start
    while repeated_end < len(input_grid) and input_grid[repeated_end] == repeated_number:
        repeated_end += 1
    
    # Create the output grid
    output_grid = input_grid[:]
    if repeated_number is not None:
        # Calculate the length of the repeated sequence
        repeated_length = repeated_end - repeated_start
        
        # Shift the repeated sequence to start after the first non-zero number
        output_grid[first_non_zero_index + 1:first_non_zero_index + 1 + repeated_length] = [repeated_number] * repeated_length
        
        # Clear the old repeated sequence
        output_grid[repeated_start:repeated_end] = [0] * repeated_length
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 2, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")