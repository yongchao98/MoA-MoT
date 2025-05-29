def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # If there's no non-zero element, return the input as is
    if first_non_zero_index is None:
        return input_grid
    
    # Create the output grid with the first non-zero element at the start
    output_grid = [input_grid[first_non_zero_index]] + input_grid[:-1]
    
    return output_grid

# Test input
input_grid = [0, 8, 0, 0, 0, 0, 0, 0, 7, 3, 0, 0, 4, 6]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")