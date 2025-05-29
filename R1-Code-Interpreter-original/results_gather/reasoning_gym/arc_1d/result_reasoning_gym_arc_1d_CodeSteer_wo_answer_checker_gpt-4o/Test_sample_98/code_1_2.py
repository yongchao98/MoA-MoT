def transform_grid(input_grid):
    output_grid = input_grid.copy()
    n = len(input_grid)
    i = 0

    while i < n:
        if input_grid[i] != 0:
            # Start of a non-zero sequence
            start = i
            while i < n and input_grid[i] == input_grid[start]:
                i += 1
            # End of the sequence
            end = i
            # Determine the length of the sequence to extend
            sequence_length = end - start
            if sequence_length == 1:
                extension_length = 3
            else:
                extension_length = max(4, sequence_length)
            for j in range(start, min(start + extension_length, n)):
                output_grid[j] = input_grid[start]
        else:
            i += 1

    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 7, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")