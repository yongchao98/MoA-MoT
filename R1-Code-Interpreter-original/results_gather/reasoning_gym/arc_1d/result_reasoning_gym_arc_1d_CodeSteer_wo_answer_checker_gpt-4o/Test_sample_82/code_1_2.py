# Test input
input_grid = [6, 8, 4, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 2]

# Find the first non-zero sequence
start_index = 0
while start_index < len(input_grid) and input_grid[start_index] == 0:
    start_index += 1

# Extract the non-zero sequence
non_zero_sequence = []
while start_index < len(input_grid) and input_grid[start_index] != 0:
    non_zero_sequence.append(input_grid[start_index])
    start_index += 1

# Create the output grid
output_grid = non_zero_sequence + input_grid[start_index:]

# Print the output grid
print("<<<", " ".join(map(str, output_grid)), ">>>")