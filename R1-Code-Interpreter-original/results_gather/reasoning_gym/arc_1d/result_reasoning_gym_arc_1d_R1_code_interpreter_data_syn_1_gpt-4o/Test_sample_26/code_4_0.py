def transform_grid(input_grid):
    # Find the first and last non-zero elements
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    last_non_zero = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1
    
    # Extract the non-zero sequence
    non_zero_sequence = input_grid[first_non_zero:last_non_zero+1]
    
    # Count zeros
    zero_count = input_grid.count(0)
    
    # Calculate the number of zeros to place at the start
    zeros_at_start = (len(input_grid) - len(non_zero_sequence)) // 2
    
    # Create the output grid
    output_grid = [0] * zeros_at_start + non_zero_sequence + [0] * (len(input_grid) - zeros_at_start - len(non_zero_sequence))
    
    return output_grid

# Test input
input_grid = [2, 2, 2, 2, 2, 2, 2, 0, 0, 2, 2, 2, 2, 2, 2, 2]
output_grid = transform_grid(input_grid)
print(output_grid)