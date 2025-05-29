def generate_output_grid(input_grid):
    output_grid = []
    for i in range(0, 11, 4):
        row = []
        for j in range(0, 11, 4):
            # Check the 3x3 block in the input grid
            block_has_one = False
            for x in range(i, min(i + 3, 11)):
                for y in range(j, min(j + 3, 11)):
                    if input_grid[x][y] == 1:
                        block_has_one = True
                        break
                if block_has_one:
                    break
            # Determine the output value for this block
            if block_has_one:
                row.append(4)
            else:
                row.append(6)
        output_grid.append(row)
    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 9, 6, 6, 6, 9, 6, 1, 6],
    [6, 1, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [6, 6, 6, 9, 6, 1, 6, 9, 1, 6, 6],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 9, 6, 6, 1, 9, 6, 6, 6],
    [6, 6, 6, 9, 6, 1, 6, 9, 6, 6, 1],
    [6, 1, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [1, 6, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [6, 6, 6, 9, 1, 6, 6, 9, 6, 1, 6],
    [1, 6, 6, 9, 6, 6, 1, 9, 6, 6, 6]
]

output_grid = generate_output_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))