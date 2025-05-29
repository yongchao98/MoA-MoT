def transform_grid(input_grid):
    # Convert the input grid to a list of integers
    grid = list(map(int, input_grid.split()))
    
    # Identify the first non-zero sequence
    non_zero_value = None
    start_index = None
    for i, value in enumerate(grid):
        if value != 0:
            non_zero_value = value
            start_index = i
            break
    
    # If no non-zero value is found, return the grid as is
    if non_zero_value is None:
        return grid
    
    # Determine the length of the non-zero sequence
    length = 0
    for i in range(start_index, len(grid)):
        if grid[i] == non_zero_value:
            length += 1
        else:
            break
    
    # Create the output grid
    output_grid = [0] * len(grid)
    
    # Place the non-zero sequence shifted by 2 positions
    for i in range(length):
        if start_index + i + 2 < len(grid):
            output_grid[start_index + i + 2] = non_zero_value
    
    # Fill the trailing positions with the non-zero value
    for i in range(len(grid) - length, len(grid)):
        output_grid[i] = non_zero_value
    
    return output_grid

# Test input
input_grid = "0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0"
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")