def expand_grid(input_grid):
    output_grid = []
    for row in input_grid:
        expanded_rows = [[element] * 3 for element in row]
        for _ in range(3):
            output_grid.append([item for sublist in expanded_rows for item in sublist])
    return output_grid

input_grid = [
    [6, 5],
    [6, 9]
]

output_grid = expand_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))