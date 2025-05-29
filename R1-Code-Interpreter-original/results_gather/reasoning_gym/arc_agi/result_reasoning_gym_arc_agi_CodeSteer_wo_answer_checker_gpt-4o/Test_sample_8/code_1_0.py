def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    for i, row in enumerate(input_grid):
        j = 0
        while j < len(row):
            if row[j] == 8:
                # Start replacing with 0 from the next element
                k = j + 1
                while k < len(row) and row[k] != 8:
                    output_grid[i][k] = 0
                    k += 1
                # Move j to the next 8 or end of row
                j = k
            else:
                j += 1

    return output_grid

# Test input grid
input_grid = [
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 8, 5, 8, 8, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 8, 5, 5, 8, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 8, 5, 5, 8, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 5, 5, 5, 5],
    [5, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 8, 5, 5, 8, 5, 5, 5, 8, 8, 8, 8, 5, 5, 5],
    [5, 8, 5, 5, 8, 5, 5, 5, 8, 5, 5, 8, 5, 5, 5],
    [5, 8, 5, 8, 8, 5, 5, 5, 8, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 8, 8, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))