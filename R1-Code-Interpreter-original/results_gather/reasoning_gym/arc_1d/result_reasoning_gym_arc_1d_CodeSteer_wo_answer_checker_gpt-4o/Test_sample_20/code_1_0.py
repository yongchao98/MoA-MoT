# Test input
input_grid = [0, 1, 2, 1, 0, 0, 1, 0, 1, 0, 2]

# Applying the derived rule
output_grid = [2, 2]  # Start with two 2s
output_grid.extend([x for x in input_grid if x == 0])  # Append zeros
output_grid.extend([x for x in input_grid if x == 1])  # Append ones

# Print the final output grid
print("<<<", " ".join(map(str, output_grid)), ">>>")