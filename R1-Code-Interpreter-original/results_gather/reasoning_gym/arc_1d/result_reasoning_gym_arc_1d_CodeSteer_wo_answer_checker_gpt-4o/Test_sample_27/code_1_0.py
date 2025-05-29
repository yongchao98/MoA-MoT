def transform_grid(input_grid):
    # Move the first element to the end of the list
    output_grid = input_grid[1:] + input_grid[:1]
    return output_grid

# Test input
input_grid = [7, 9, 6, 4, 0, 0, 0, 8, 7, 1, 3, 5]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")