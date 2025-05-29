def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    for i, row in enumerate(input_grid):
        if 2 in row:
            # Find the start and end indices of the sequence of '2's
            start = row.index(2)
            end = len(row) - row[::-1].index(2) - 1

            # Determine the replacement based on the context
            if start > 0 and end < len(row) - 1:
                # If '2's are surrounded by the same number, use that number
                if row[start - 1] == row[end + 1]:
                    replacement = row[start - 1]
                else:
                    # Otherwise, choose the most frequent number in the surrounding context
                    context = row[max(0, start - 2):min(len(row), end + 3)]
                    replacement = max(set(context), key=context.count)
            elif start > 0:
                replacement = row[start - 1]
            else:
                replacement = row[end + 1]

            # Replace the sequence of '2's with the determined replacement
            for j in range(start, end + 1):
                output_grid[i][j] = replacement

    return output_grid

# Test input grid
input_grid = [
    [3, 3, 1, 8, 0, 1, 5, 4, 4, 5, 1, 0, 8, 1, 3, 3],
    [3, 8, 8, 1, 1, 1, 4, 4, 4, 4, 1, 1, 1, 8, 8, 3],
    [1, 8, 4, 0, 5, 4, 4, 2, 2, 2, 4, 5, 0, 4, 8, 1],
    [8, 1, 0, 8, 4, 4, 8, 2, 2, 2, 4, 4, 8, 0, 1, 8],
    [0, 1, 2, 2, 2, 2, 4, 2, 2, 2, 3, 1, 4, 5, 1, 0],
    [1, 1, 2, 2, 2, 2, 4, 7, 7, 4, 3, 3, 4, 4, 1, 1],
    [5, 4, 2, 2, 2, 2, 7, 7, 7, 7, 4, 4, 8, 4, 4, 5],
    [4, 4, 2, 2, 2, 2, 7, 3, 3, 7, 7, 4, 8, 8, 4, 4],
    [4, 4, 8, 8, 4, 7, 7, 3, 3, 7, 7, 4, 8, 8, 4, 4],
    [5, 4, 4, 8, 4, 4, 7, 7, 7, 7, 4, 4, 8, 4, 4, 5],
    [1, 1, 4, 4, 3, 3, 4, 7, 7, 4, 3, 3, 4, 4, 1, 1],
    [0, 1, 5, 4, 1, 3, 4, 4, 4, 4, 3, 1, 4, 5, 1, 0],
    [8, 1, 0, 8, 4, 4, 8, 8, 8, 8, 4, 4, 8, 0, 1, 8],
    [1, 8, 4, 0, 5, 4, 4, 8, 8, 4, 4, 5, 0, 4, 8, 1],
    [3, 8, 8, 1, 1, 1, 4, 4, 4, 4, 1, 1, 1, 8, 8, 3],
    [3, 3, 1, 8, 0, 1, 5, 4, 4, 5, 1, 0, 8, 1, 3, 3]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))