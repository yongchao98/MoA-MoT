def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is not None:
        # Find the sequence of the next number
        sequence_start = first_non_zero_index + 1
        sequence_number = input_grid[sequence_start]
        
        # Check if there is a sequence to shift
        if sequence_number != 0:
            # Shift the sequence left by one position
            for i in range(sequence_start, len(input_grid) - 1):
                input_grid[i] = input_grid[i + 1]
            # Replace the last number in the sequence with 0
            input_grid[-1] = 0
    
    return input_grid

# Test input
test_input = [9, 9, 9, 9, 9, 9, 9, 9, 0, 0, 0, 2]
# Transform the test input
output_grid = transform_grid(test_input)
print(output_grid)