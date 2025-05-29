def transform_grid(input_grid):
    # Find the first non-zero element from the left
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    # Find the first zero element from the right
    first_zero_index_from_right = next((i for i, x in enumerate(reversed(input_grid)) if x == 0), None)
    
    if first_non_zero_index is not None and first_zero_index_from_right is not None:
        # Calculate the actual index from the right
        first_zero_index_from_right = len(input_grid) - 1 - first_zero_index_from_right
        
        # If zeros are at the start, shift them to the right
        if first_non_zero_index > 0:
            output_grid = input_grid[first_non_zero_index:] + input_grid[:first_non_zero_index]
        # If zeros are at the end, shift them to the left
        elif first_zero_index_from_right < len(input_grid) - 1:
            output_grid = input_grid[first_zero_index_from_right + 1:] + input_grid[:first_zero_index_from_right + 1]
        else:
            output_grid = input_grid
    else:
        output_grid = input_grid
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)