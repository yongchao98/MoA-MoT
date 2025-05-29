def extract_7x7_block(input_grid):
    # Based on the examples, the starting point for the 7x7 block seems to be around the middle of the grid.
    # Let's assume the starting point is (4, 4) for this test input.
    start_row, start_col = 4, 4
    output_grid = []
    for i in range(7):
        row = input_grid[start_row + i][start_col:start_col + 7]
        output_grid.append(row)
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 1, 9, 9, 9, 0, 0, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 1, 9, 1, 9, 9, 0, 9, 9, 9],
    [9, 9, 9, 8, 9, 8, 9, 9, 1, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 8, 9, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 7, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 5, 6, 6, 6, 8, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 5, 9, 5, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 5, 5, 5, 9, 9, 9, 9],
    [6, 7, 6, 6, 6, 0, 6, 9, 9, 9, 5, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9]
]

output_grid = extract_7x7_block(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))