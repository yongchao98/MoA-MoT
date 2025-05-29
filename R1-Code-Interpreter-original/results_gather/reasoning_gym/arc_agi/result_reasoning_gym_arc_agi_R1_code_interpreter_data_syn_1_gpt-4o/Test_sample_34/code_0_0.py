def transform_grid(input_grid):
    output_grid = []
    # First two rows remain unchanged
    output_grid.append(input_grid[0])
    output_grid.append(input_grid[1])
    # Third row becomes all 2s
    output_grid.append([2] * len(input_grid[2]))
    # Fourth row remains unchanged
    output_grid.append(input_grid[3])
    # Fifth row alternates 9 and 1
    fifth_row = input_grid[4]
    new_fifth_row = []
    toggle = True
    for value in fifth_row:
        if value == 9:
            new_fifth_row.append(1 if toggle else 9)
            toggle = not toggle
        else:
            new_fifth_row.append(value)
    output_grid.append(new_fifth_row)
    return output_grid

# Test input grid
input_grid = [
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
    [1, 2, 1, 2, 2],
    [9, 2, 9, 2, 2],
    [9, 9, 9, 9, 9]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))