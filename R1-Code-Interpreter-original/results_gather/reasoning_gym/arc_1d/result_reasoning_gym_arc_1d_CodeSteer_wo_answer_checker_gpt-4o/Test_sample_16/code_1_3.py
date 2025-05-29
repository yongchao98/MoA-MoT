def transform_grid(input_grid):
    # Shift all elements to the right by one position
    output_grid = [0] + input_grid[:-1]
    
    return output_grid

# Test input
input_grid = [0, 8, 0, 0, 0, 0, 0, 0, 7, 3, 0, 0, 4, 6]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")