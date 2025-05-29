def transform_grid(input_grid):
    n = len(input_grid)
    output_grid = [0] * n
    i = 0

    while i < n:
        if input_grid[i] != 0:
            # Start of a sequence
            start = i
            while i < n and input_grid[i] == input_grid[start]:
                i += 1
            end = i

            # Check if the sequence should be retained
            if i < n and input_grid[i] == 0:
                # Retain the sequence
                for j in range(start, end):
                    output_grid[j] = input_grid[start]
        else:
            i += 1

    return output_grid

# Test input
input_grid = [0, 5, 0, 5, 0, 0, 0, 0, 0, 0, 5, 0, 0, 5, 5, 5, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")