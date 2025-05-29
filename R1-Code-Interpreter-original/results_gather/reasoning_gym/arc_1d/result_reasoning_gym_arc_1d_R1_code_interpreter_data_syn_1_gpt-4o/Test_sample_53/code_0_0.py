def transform_grid(input_grid):
    # Count the number of zeros
    zero_count = input_grid.count(0)
    # Count the number of non-zeros
    non_zero_count = len(input_grid) - zero_count
    # Create the output grid with zeros on the left and non-zeros on the right
    output_grid = [0] * zero_count + [7] * non_zero_count
    return output_grid

# Test input
input_grid = [0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)