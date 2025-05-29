def transform_grid(input_grid):
    output_grid = input_grid[:]
    n = len(input_grid)
    
    for i in range(n):
        if input_grid[i] != 0:
            j = i + 1
            while j < n and input_grid[j] == 0:
                output_grid[j] = input_grid[i]
                j += 1
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 0, 9]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")