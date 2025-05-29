def transform_row(row):
    # Find the first occurrence of 7
    if 7 in row:
        index_of_7 = row.index(7)
        # Find the number before the first 7
        if index_of_7 > 0:
            number_before_7 = row[index_of_7 - 1]
            # Find the first occurrence of 8 after the 7
            for i in range(index_of_7 + 1, len(row)):
                if row[i] == 8:
                    row[i] = number_before_7
                    break
    return row

# Test input grid
input_grid = [
    [9, 9, 8, 2, 9, 9, 9],
    [6, 9, 9, 2, 9, 9, 9],
    [9, 9, 9, 2, 9, 9, 9]
]

# Transform each row of the input grid
output_grid = [transform_row(row) for row in input_grid]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))