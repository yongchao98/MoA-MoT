def find_significant_starting_point(input_grid):
    # Find a strategic starting point by analyzing clusters of non-9 values
    for i in range(len(input_grid) - 6):
        for j in range(len(input_grid[i]) - 6):
            # Check if the 7x7 block starting at (i, j) contains significant non-9 values
            block = [input_grid[x][j:j+7] for x in range(i, i+7)]
            if any(9 not in row for row in block):
                return i, j
    return 0, 0  # Default to (0, 0) if no significant block is found

def extract_7x7_block(input_grid):
    start_row, start_col = find_significant_starting_point(input_grid)
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