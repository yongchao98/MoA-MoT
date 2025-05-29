def transform_grid(input_grid):
    return [5 if x != 0 else 0 for x in input_grid]

# Test input
input_grid = [0, 0, 0, 0, 8, 0, 7, 9, 7, 0]
output_grid = transform_grid(input_grid)

# Print the output grid
print("<<<", " ".join(map(str, output_grid)), ">>>")