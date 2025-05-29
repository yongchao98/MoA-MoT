def transform_grid(input_grid):
    output_grid = input_grid[:]
    n = len(input_grid)
    
    # Find the first non-zero value and its index
    non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    
    if non_zero_index is not None:
        value = input_grid[non_zero_index]
        
        if value == 2:
            # Fill every third position starting from the first zero before the non-zero value
            start = non_zero_index
            while start > 0 and input_grid[start - 1] == 0:
                start -= 1
            for i in range(start, n, 3):
                if input_grid[i] == 0:
                    output_grid[i] = value
                else:
                    break
        
        elif value in [5, 8]:
            # Fill every position up to the first zero after the non-zero value
            for i in range(non_zero_index, n):
                if input_grid[i] == 0:
                    output_grid[i] = value
                else:
                    break
    
    return output_grid

# Test input
input_grid = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")