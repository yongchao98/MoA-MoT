def transform_grid(input_grid):
    output_grid = input_grid[:]
    i = 0
    while i < len(input_grid):
        if input_grid[i] == 2:
            # Count the length of the contiguous block of '2's
            start = i
            while i < len(input_grid) and input_grid[i] == 2:
                i += 1
            end = i
            # Transform the first two '2's in the block to '8's
            for j in range(start, min(start + 2, end)):
                output_grid[j] = 8
        else:
            i += 1
    return output_grid

# Test input
input_grid = [0, 0, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 2, 2, 2]
output_grid = transform_grid(input_grid)
print(output_grid)