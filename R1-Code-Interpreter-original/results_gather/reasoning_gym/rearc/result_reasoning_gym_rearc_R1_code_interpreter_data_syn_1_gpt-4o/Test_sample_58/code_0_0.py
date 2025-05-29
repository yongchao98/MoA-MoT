def transform_grid(input_grid):
    # Find the row with a unique number
    for i, row in enumerate(input_grid):
        unique_numbers = set(row)
        if len(unique_numbers) < len(row):  # There is a repeated number
            unique_number = min(unique_numbers, key=row.count)
            # Extend this number to the right
            for j in range(len(row)):
                if row[j] == unique_number:
                    for k in range(j, len(row)):
                        if row[k] != unique_number:
                            break
                        row[k] = unique_number
            break
    return input_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1],
    [3, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 2, 1]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))