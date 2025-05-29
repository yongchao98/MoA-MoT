def transform_row(row):
    # Find the majority number
    majority = max(set(row), key=row.count)
    # Find the first number that is different from the majority
    for num in row:
        if num != majority:
            different_num = num
            break
    # Replace the first occurrence of the majority number after the different number
    output_row = row[:]
    different_found = False
    for i in range(len(row)):
        if row[i] == different_num:
            different_found = True
        if different_found and row[i] == majority:
            output_row[i] = different_num
            break
    return output_row

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