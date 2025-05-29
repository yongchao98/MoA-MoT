def transform_grid(input_grid):
    # Find the index of the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If there is no non-zero element, return the input grid as is
    if first_non_zero_index is None:
        return input_grid
    
    # Create the output grid by shifting the first non-zero element and subsequent elements
    output_grid = input_grid[:first_non_zero_index] + [0] + input_grid[first_non_zero_index:-1]
    
    return output_grid

# Test input grid
input_grid = [0, 9, 2, 1, 2, 8, 6, 6, 9, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")