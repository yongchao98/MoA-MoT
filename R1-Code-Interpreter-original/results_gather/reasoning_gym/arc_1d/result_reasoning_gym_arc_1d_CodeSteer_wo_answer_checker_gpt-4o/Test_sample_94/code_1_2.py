def transform_grid(input_grid):
    # Count leading zeros
    leading_zeros = 0
    for num in input_grid:
        if num == 0:
            leading_zeros += 1
        else:
            break
    
    # Shift non-zero numbers to the right by the number of leading zeros
    non_zero_part = input_grid[leading_zeros:]
    transformed_grid = [0] * leading_zeros + non_zero_part + [0] * leading_zeros
    
    # Ensure the output grid has the same length as the input grid
    transformed_grid = transformed_grid[:len(input_grid)]
    
    return transformed_grid

# Test input
input_grid = [0, 0, 8, 3, 3, 1, 9, 2, 5, 7, 5, 2, 9, 2, 8, 3, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")