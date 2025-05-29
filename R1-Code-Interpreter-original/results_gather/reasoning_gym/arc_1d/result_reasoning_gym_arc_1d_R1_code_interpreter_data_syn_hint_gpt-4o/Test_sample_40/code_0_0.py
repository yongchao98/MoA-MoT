def transform_grid(input_grid):
    # Find the first and last non-zero elements
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    last_non_zero = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero:last_non_zero + 1]
    
    # Calculate the shift amount (3 positions to the left)
    shift_amount = 3
    
    # Create the output grid
    output_grid = [0] * len(input_grid)
    output_grid[shift_amount:shift_amount + len(non_zero_sequence)] = non_zero_sequence
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)