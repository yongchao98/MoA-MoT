# Test input
input_grid = [6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 6]

# Find the first occurrence of 0
first_zero_index = input_grid.index(0)

# Construct the output grid by shifting the first '0' one position to the right
output_grid = (
    input_grid[:first_zero_index] +  # Elements before the first '0'
    [input_grid[first_zero_index + 1]] +  # The element after the first '0'
    [0] +  # The shifted '0'
    input_grid[first_zero_index + 2:]  # Elements after the shifted '0'
)

# Print the output grid
print("<<<", " ".join(map(str, output_grid)), ">>>")