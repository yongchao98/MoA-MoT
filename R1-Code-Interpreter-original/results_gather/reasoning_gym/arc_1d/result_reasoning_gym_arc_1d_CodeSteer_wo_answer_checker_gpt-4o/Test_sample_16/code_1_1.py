def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    # Find the first zero position
    first_zero_index = next((i for i, x in enumerate(input_grid) if x == 0), None)
    
    # If there's no non-zero element or no zero, return the input as is
    if first_non_zero_index is None or first_zero_index is None:
        return input_grid
    
    # Create the output grid
    output_grid = input_grid[:]
    
    # Move the first non-zero element to the first zero position
    output_grid[first_zero_index] = input_grid[first_non_zero_index]
    
    # Shift elements before the first zero to the right
    for i in range(first_zero_index):
        output_grid[i] = input_grid[i + 1]
    
    return output_grid

# Test input
input_grid = [0, 8, 0, 0, 0, 0, 0, 0, 7, 3, 0, 0, 4, 6]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")