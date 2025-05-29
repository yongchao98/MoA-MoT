def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        unique_elements = set(row)
        if len(unique_elements) == 2:
            # Find the element that is different
            for element in unique_elements:
                if row.count(element) == 1:
                    output_grid.append([element] * len(row))
                    break
        else:
            output_grid.append(row)
    return output_grid

# Test input grid
input_grid = [
    [0, 3, 3, 3, 3],
    [3, 3, 3, 3, 3],
    [3, 3, 0, 3, 3],
    [3, 3, 3, 3, 3],
    [3, 3, 3, 3, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))