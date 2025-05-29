def find_output_grid(input_grid):
    for row in input_grid:
        if 1 in row:
            index_of_1 = row.index(1)
            if index_of_1 >= 2:
                # Check if the numbers before '1' are identical
                if row[index_of_1 - 2] == row[index_of_1 - 1]:
                    num = row[index_of_1 - 1]
                    # Return the 2x2 grid
                    return [[num, num], [num, num]]

input_grid = [
    [8, 8, 8, 8, 8, 1, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 1, 8, 8, 8, 8, 8],
    [8, 9, 9, 9, 9, 1, 8, 8, 8, 9, 9],
    [8, 8, 8, 8, 8, 1, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 1, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 1, 8, 8, 8, 8, 8],
    [8, 4, 4, 4, 4, 1, 8, 8, 4, 4, 4],
    [8, 8, 8, 8, 8, 1, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 1, 8, 8, 8, 8, 8]
]

output_grid = find_output_grid(input_grid)
print(output_grid)