def transform_grid(input_grid):
    n = len(input_grid)
    output_grid = [0] * n
    
    # Use a sliding window to find the sequence of non-zero numbers
    start = 0
    while start < n and input_grid[start] == 0:
        start += 1
    
    if start == n:
        return input_grid  # No non-zero numbers found
    
    end = start
    while end < n and input_grid[end] != 0:
        end += 1
    
    # Copy the sequence to the output grid
    output_grid[start:end] = input_grid[start:end]
    
    # If the sequence is not at the end, shift it to the left
    if end < n:
        shift_length = end - start
        output_grid[start:start + shift_length] = input_grid[start:end]
    
    return output_grid

# Test input
input_grid = [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")