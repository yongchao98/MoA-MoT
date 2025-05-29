def transform_grid(input_grid):
    # First 10 numbers are replaced by zeros
    output_grid = [0] * 10
    
    # Next 4 numbers are moved to the start of the output after the zeros
    output_grid.extend(input_grid[10:14])
    
    # Remaining numbers are shifted to the left by 10 positions
    output_grid.extend(input_grid[14:])
    
    # First 4 numbers of the input are appended to the end of the output
    output_grid.extend(input_grid[:4])
    
    return output_grid

# Test input
input_grid = [5, 5, 9, 7, 4, 1, 5, 6, 9, 2, 7, 8, 1, 7, 1, 6, 8, 9, 9, 8, 8, 5, 4, 7, 4, 9, 6, 7, 7]
output_grid = transform_grid(input_grid)
print(output_grid)