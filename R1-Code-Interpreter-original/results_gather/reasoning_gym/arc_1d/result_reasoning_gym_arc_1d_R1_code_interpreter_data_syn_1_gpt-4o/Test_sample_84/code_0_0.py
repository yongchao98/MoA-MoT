def transform_grid(input_grid):
    output_grid = []
    i = 0
    while i < len(input_grid):
        output_grid.append(input_grid[i])
        if input_grid[i] != 0:
            # Check if it's the end of a sequence
            if i == len(input_grid) - 1 or input_grid[i + 1] == 0:
                output_grid.append(3)
        i += 1
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 5, 5]
output_grid = transform_grid(input_grid)
print(output_grid)