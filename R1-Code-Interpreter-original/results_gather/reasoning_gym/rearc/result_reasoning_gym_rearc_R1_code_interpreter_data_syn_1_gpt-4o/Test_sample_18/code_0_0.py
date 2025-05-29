def find_repeated_row(input_grid):
    # Convert each row to a tuple to make it hashable
    row_counts = {}
    for row in input_grid:
        row_tuple = tuple(row)
        if row_tuple in row_counts:
            row_counts[row_tuple] += 1
        else:
            row_counts[row_tuple] = 1
    
    # Find the row with the minimum count (most unique)
    min_count = min(row_counts.values())
    for row, count in row_counts.items():
        if count == min_count:
            return list(row)

def transform_grid(input_grid):
    repeated_row = find_repeated_row(input_grid)
    output_grid = []
    for i, row in enumerate(input_grid):
        if i % 2 == 0:
            output_grid.append(repeated_row)
        else:
            output_grid.append(row)
    return output_grid

# Test input grid
input_grid = [
    [7, 0, 9, 0, 9],
    [0, 0, 0, 0, 0],
    [9, 0, 9, 0, 9],
    [0, 0, 0, 0, 0],
    [7, 0, 9, 0, 9],
    [0, 0, 0, 0, 0],
    [9, 0, 9, 0, 9]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))