def apply_shift_rule(input_grid):
    # Count leading zeros
    leading_zeros = 0
    for num in input_grid:
        if num == 0:
            leading_zeros += 1
        else:
            break
    
    # Shift non-zero elements to the right
    non_zero_elements = [num for num in input_grid if num != 0]
    output_grid = [0] * leading_zeros + non_zero_elements
    
    # Ensure the output grid has the same length as the input grid
    output_grid = output_grid[:len(input_grid)]
    
    return output_grid

# Test input grid
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 8, 9, 4, 2, 5, 4, 7, 7, 2, 9, 7, 0, 0, 0]
output_grid = apply_shift_rule(input_grid)
print(output_grid)