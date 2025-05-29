def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if first_non_zero_index is not None:
        # Create the output grid by shifting elements
        output_grid = [0] * len(input_grid)
        output_grid[0] = input_grid[first_non_zero_index]
        output_grid[1:] = input_grid[:first_non_zero_index] + input_grid[first_non_zero_index+1:]
    else:
        # If no non-zero element is found, return the input grid as is
        output_grid = input_grid
    
    return output_grid

# Test input
input_grid = [0, 8, 0, 0, 0, 0, 0, 0, 7, 3, 0, 0, 4, 6]
output_grid = transform_grid(input_grid)
print(output_grid)