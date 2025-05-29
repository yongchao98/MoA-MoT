def transform_grid(input_grid):
    # Function to determine the replacement for a sequence of '2's
    def determine_replacement(row, start, end):
        # Analyze the context around the sequence of '2's
        context = row[max(0, start - 2):min(len(row), end + 3)]
        # Use a dictionary to count occurrences of each number in the context
        count = {}
        for num in context:
            if num != 2:  # Ignore '2's in the context
                count[num] = count.get(num, 0) + 1
        # Find the most frequent number in the context
        replacement = max(count, key=count.get)
        return replacement

    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    for i, row in enumerate(input_grid):
        if 2 in row:
            # Find the start and end indices of the sequence of '2's
            start = row.index(2)
            end = len(row) - row[::-1].index(2) - 1

            # Determine the replacement based on the context
            replacement = determine_replacement(row, start, end)

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