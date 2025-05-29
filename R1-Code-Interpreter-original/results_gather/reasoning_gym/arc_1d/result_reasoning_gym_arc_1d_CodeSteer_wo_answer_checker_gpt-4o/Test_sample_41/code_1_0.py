def transform_grid(input_grid):
    # Find the first occurrence of zero
    first_zero_index = input_grid.index(0)
    # Remove the first zero and append a zero at the end
    output_grid = input_grid[:first_zero_index] + input_grid[first_zero_index+1:] + [0]
    return output_grid

# Test input
input_grid = [3, 8, 4, 0, 0, 0, 7, 7, 1, 0, 5, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<", " ".join(map(str, output_grid)), ">>>")