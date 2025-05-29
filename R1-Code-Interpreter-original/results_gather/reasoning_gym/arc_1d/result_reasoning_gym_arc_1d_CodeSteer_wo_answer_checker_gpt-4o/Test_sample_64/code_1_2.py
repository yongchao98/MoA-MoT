def transform_grid(input_grid):
    # Find the first non-zero number and its index
    first_non_zero_index = None
    for i, num in enumerate(input_grid):
        if num != 0:
            first_non_zero_index = i
            break
    
    # If no non-zero number is found, return the input grid as is
    if first_non_zero_index is None:
        return input_grid
    
    # Determine the number that follows the first non-zero number
    replacement_number = input_grid[first_non_zero_index + 1]
    
    # Replace the first non-zero number and all subsequent numbers
    output_grid = input_grid[:first_non_zero_index] + [replacement_number] * (len(input_grid) - first_non_zero_index)
    
    return output_grid

# Test input
input_grid = [2, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")