def transform_grid(input_grid):
    output_grid = input_grid[:]
    n = len(input_grid)
    
    # Find the first non-zero value and its index
    non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if non_zero_index is not None:
        value = input_grid[non_zero_index]
        
        # Determine the fill range
        left = non_zero_index
        right = non_zero_index
        
        # Expand to the left
        while left > 0 and input_grid[left - 1] == 0:
            left -= 1
        
        # Expand to the right
        while right < n - 1 and input_grid[right + 1] == 0:
            right += 1
        
        # Fill the determined range with the non-zero value
        for i in range(left, right + 1):
            output_grid[i] = value
    
    return output_grid

# Test input
input_grid = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")