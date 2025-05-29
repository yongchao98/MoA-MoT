def transform_grid(input_grid):
    # Identify the non-zero sequence
    non_zero_sequence = [x for x in input_grid if x != 0]
    
    # Determine the shift amount (4 positions to the right)
    shift_amount = 4
    
    # Create the output grid with zeros
    output_grid = [0] * len(input_grid)
    
    # Place the non-zero sequence in the new position
    for i, value in enumerate(non_zero_sequence):
        output_grid[i + shift_amount] = value
    
    return output_grid

# Test input
input_grid = [9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9]
output_grid = transform_grid(input_grid)
print(output_grid)