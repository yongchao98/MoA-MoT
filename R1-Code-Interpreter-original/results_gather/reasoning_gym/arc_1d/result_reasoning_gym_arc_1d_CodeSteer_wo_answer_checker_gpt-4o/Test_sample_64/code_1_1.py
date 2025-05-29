def transform_grid(input_grid):
    # Find the first non-zero number
    for i, num in enumerate(input_grid):
        if num != 0:
            # Replace from this index onwards with the number that follows the first non-zero number
            return input_grid[:i] + [input_grid[i+1]] * (len(input_grid) - i)
    return input_grid

# Test input
input_grid = [2, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")