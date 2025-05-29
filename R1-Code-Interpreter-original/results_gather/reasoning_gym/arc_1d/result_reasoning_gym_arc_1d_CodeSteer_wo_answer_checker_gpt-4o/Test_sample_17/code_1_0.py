def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if first_non_zero is None:
        return input_grid  # If all are zeros, return as is

    # Find the last non-zero element
    last_non_zero = len(input_grid) - 1 - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)

    # Shift the non-zero sequence to the left by one position
    output_grid = input_grid[:first_non_zero] + input_grid[first_non_zero+1:last_non_zero+1] + [0] + input_grid[last_non_zero+1:]

    return output_grid

# Test input
input_grid = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")