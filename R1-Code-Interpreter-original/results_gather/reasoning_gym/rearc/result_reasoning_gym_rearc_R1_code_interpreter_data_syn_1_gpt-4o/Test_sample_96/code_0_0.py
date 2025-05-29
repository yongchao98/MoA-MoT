def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        # Count the frequency of each number in the row, excluding zeros
        freq = {}
        for num in row:
            if num != 0:
                if num in freq:
                    freq[num] += 1
                else:
                    freq[num] = 1
        # Find the most frequent number
        most_frequent = max(freq, key=freq.get)
        # Replace zeros with the most frequent number
        new_row = [num if num != 0 else most_frequent for num in row]
        output_grid.append(new_row)
    return output_grid

# Test input grid
input_grid = [
    [4, 0, 4, 0, 4, 0, 4, 0, 4],
    [0, 4, 0, 4, 0, 4, 0, 4, 0],
    [4, 0, 4, 0, 4, 0, 4, 0, 4],
    [0, 4, 0, 4, 0, 4, 0, 4, 0],
    [4, 0, 4, 0, 4, 0, 4, 0, 4],
    [6, 6, 0, 4, 0, 4, 0, 4, 0],
    [6, 6, 4, 0, 4, 0, 4, 0, 4],
    [0, 4, 0, 4, 0, 4, 0, 4, 0],
    [4, 0, 4, 0, 4, 0, 4, 0, 4]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))