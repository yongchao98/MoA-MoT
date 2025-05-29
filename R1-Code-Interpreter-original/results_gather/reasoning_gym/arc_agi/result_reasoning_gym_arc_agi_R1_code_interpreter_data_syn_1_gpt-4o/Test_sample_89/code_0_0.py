# Test input grid
input_grid = [
    [4, 4, 5, 5],
    [5, 4, 4, 4],
    [5, 4, 4, 5],
    [5, 4, 4, 4],
    [4, 5, 4, 4],
    [4, 5, 4, 5],
    [5, 4, 5, 5],
    [4, 4, 4, 5],
    [4, 5, 4, 5],
    [5, 5, 5, 4],
    [4, 5, 4, 5],
    [4, 4, 4, 5],
    [4, 5, 5, 5],
    [4, 4, 5, 5],
    [4, 4, 5, 5],
    [5, 5, 5, 4],
    [5, 4, 4, 5],
    [5, 4, 4, 4],
    [5, 4, 4, 5],
    [5, 5, 4, 4],
    [4, 5, 5, 5],
    [4, 5, 5, 5],
    [5, 5, 5, 5],
    [5, 5, 4, 5]
]

# Applying the transformation rule
output_grid = []
for i, row in enumerate(input_grid):
    if 12 <= i <= 16:  # Based on the pattern observed in examples
        # Change the first two elements to 7
        new_row = [7, 7] + row[2:]
    else:
        new_row = row
    output_grid.append(new_row)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))