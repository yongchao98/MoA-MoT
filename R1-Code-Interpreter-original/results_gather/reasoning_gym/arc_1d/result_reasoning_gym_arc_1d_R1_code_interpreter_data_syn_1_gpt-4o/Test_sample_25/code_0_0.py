def transform_grid(input_grid):
    # Count leading zeros
    leading_zeros = 0
    for num in input_grid:
        if num == 0:
            leading_zeros += 1
        else:
            break
    
    # Shift the non-zero sequence to the right by the number of leading zeros
    non_zero_sequence = input_grid[leading_zeros:]
    output_grid = [0] * leading_zeros + [0] * leading_zeros + non_zero_sequence[:-leading_zeros]
    
    return output_grid

# Test input grid
input_grid = [0, 0, 0, 1, 8, 8, 1, 7, 1, 6, 8, 9, 1, 8, 9, 1, 5, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)