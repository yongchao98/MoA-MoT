def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
    
    # Shift non-zero elements to the left
    non_zero_part = input_grid[first_non_zero_index:]
    
    # Count leading zeros
    leading_zeros_count = first_non_zero_index
    
    # Append leading zeros to the end
    output_grid = non_zero_part + [0] * leading_zeros_count
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 8, 9, 2, 8, 8, 6, 3, 9, 8, 2, 6, 4, 7, 7]
output_grid = transform_grid(input_grid)
print(output_grid)