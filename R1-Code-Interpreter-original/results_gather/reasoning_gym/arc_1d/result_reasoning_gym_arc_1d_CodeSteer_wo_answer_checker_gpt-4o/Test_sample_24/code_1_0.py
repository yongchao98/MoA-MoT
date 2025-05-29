def transform_grid(input_grid):
    # Count leading zeros
    leading_zeros = 0
    for num in input_grid:
        if num == 0:
            leading_zeros += 1
        else:
            break
    
    # Find the non-zero sequence
    non_zero_sequence = input_grid[leading_zeros:]
    non_zero_sequence = [num for num in non_zero_sequence if num != 0]
    
    # Create the output grid
    output_grid = [0] * leading_zeros + non_zero_sequence + [0] * (len(input_grid) - leading_zeros - len(non_zero_sequence))
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")