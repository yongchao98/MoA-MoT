def transform_grid(input_grid):
    n = len(input_grid)
    start, end = None, None
    
    # Identify the start and end of the non-zero sequence
    for i in range(n):
        if input_grid[i] != 0:
            if start is None:
                start = i
            end = i
    
    # Create the output grid with zeros
    output_grid = [0] * n
    
    # If a non-zero sequence is found, shift it two positions to the left
    if start is not None and end is not None:
        length = end - start + 1
        output_grid[2:2 + length] = input_grid[start:end + 1]
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")