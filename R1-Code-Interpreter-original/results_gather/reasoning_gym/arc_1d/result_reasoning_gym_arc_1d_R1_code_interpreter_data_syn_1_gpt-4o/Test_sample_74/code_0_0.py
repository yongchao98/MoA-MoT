# Test input grid
input_grid = [0, 0, 0, 0, 8, 0, 7, 9, 7, 0]

# Apply the rule: transform non-zero numbers to 5, keep 0 unchanged
output_grid = [5 if x != 0 else 0 for x in input_grid]

# Print the output grid
print(output_grid)