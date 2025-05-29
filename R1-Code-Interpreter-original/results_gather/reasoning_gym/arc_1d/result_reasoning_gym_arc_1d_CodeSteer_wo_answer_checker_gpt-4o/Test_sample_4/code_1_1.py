def transform_grid(input_grid):
    # Find the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is None:
        return input_grid  # If all are zeros, return as is
    
    # Extract the sequence after the first non-zero number
    sequence = input_grid[first_non_zero_index:]
    
    # Check if there are zeros in the sequence
    if 0 in sequence:
        # Remove zeros and append them to the end
        non_zero_sequence = [x for x in sequence if x != 0]
        zero_count = len(sequence) - len(non_zero_sequence)
        transformed_sequence = non_zero_sequence + [0] * zero_count
        return input_grid[:first_non_zero_index] + transformed_sequence
    else:
        return input_grid

# Test input
input_grid = [8, 8, 8, 8, 0, 0, 1, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")