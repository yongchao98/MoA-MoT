def transform_grid(input_grid):
    # Find the row with a unique number
    for i, row in enumerate(input_grid):
        # Count occurrences of each number in the row
        number_counts = {num: row.count(num) for num in set(row)}
        # Find the number with the minimum count
        unique_number = min(number_counts, key=number_counts.get)
        # Check if this number is indeed unique (appears less frequently)
        if number_counts[unique_number] < len(row):
            # Extend this number to the right
            for j in range(len(row)):
                if row[j] == unique_number:
                    for k in range(j, len(row)):
                        if row[k] != unique_number:
                            row[k] = unique_number
                        else:
                            break
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