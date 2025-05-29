def transform_grid(input_grid):
    output_grid = []
    n = len(input_grid)
    
    i = 0
    while i < n:
        output_grid.append(input_grid[i])
        
        # Check if the current number is non-zero and followed by a zero or end of grid
        if input_grid[i] != 0:
            if i + 1 < n and input_grid[i + 1] == 0:
                output_grid.append(3)
            elif i + 1 == n:
                output_grid.append(3)
        
        i += 1
    
    return output_grid

# Test input
input_grid = [4, 0, 0, 5, 5, 0, 0, 0, 9, 9, 9, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")