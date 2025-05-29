def transform_grid(input_grid):
    n = len(input_grid)
    output_grid = [0] * n
    i = 0

    while i < n:
        if input_grid[i] != 0:
            # Start of a non-zero sequence
            start = i
            while i < n and input_grid[i] != 0:
                i += 1
            end = i

            # Extract the sequence
            sequence = input_grid[start:end]

            # Determine the new position
            if start == 0:
                # If the sequence starts at the beginning, move it to the start of the output
                output_grid[:len(sequence)] = sequence
            else:
                # Otherwise, shift the sequence to the right by four positions
                new_start = min(start + 4, n - len(sequence))
                output_grid[new_start:new_start + len(sequence)] = sequence
        else:
            i += 1

    return output_grid

# Test input
input_grid = [9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")