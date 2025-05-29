def transform_grid(input_grid):
    # Identify the non-zero sequence
    non_zero_sequence = [x for x in input_grid if x != 0]
    
    # Determine the shift amount based on the position of the first non-zero number
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If the first non-zero number is in the first half, shift right by 4
    # If in the second half, shift left by 3
    if first_non_zero_index < len(input_grid) // 2:
        shift_amount = 4
    else:
        shift_amount = -3
    
    # Create the output grid with zeros
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence in the new position
    for i, value in enumerate(non_zero_sequence):
        new_index = i + shift_amount
        if 0 <= new_index < len(output_grid):
            output_grid[new_index] = value
    
    return output_grid

# Test input
input_grid = [9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9]
output_grid = transform_grid(input_grid)
print(output_grid)