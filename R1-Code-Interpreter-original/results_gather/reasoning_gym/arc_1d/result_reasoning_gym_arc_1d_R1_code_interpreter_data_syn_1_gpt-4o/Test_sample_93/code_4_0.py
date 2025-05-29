def transform_grid(input_grid):
    # Find the last non-zero sequence
    non_zero_sequence = []
    for num in reversed(input_grid):
        if num != 0:
            non_zero_sequence.append(num)
        elif non_zero_sequence:
            break
    
    # Reverse the non-zero sequence to maintain the original order
    non_zero_sequence.reverse()
    
    # Create the output grid
    output_grid = non_zero_sequence + input_grid[:-len(non_zero_sequence)]
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 2, 7, 0]
output_grid = transform_grid(input_grid)
print(output_grid)