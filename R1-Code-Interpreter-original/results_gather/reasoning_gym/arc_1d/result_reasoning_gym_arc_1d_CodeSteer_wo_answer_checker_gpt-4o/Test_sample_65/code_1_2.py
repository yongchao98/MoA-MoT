def transform_grid(input_grid):
    # Find the first occurrence of 1
    try:
        first_one_index = input_grid.index(1)
    except ValueError:
        # If there are no 1s, shift the sequence of 0s
        return input_grid[1:] + [input_grid[0]]

    # Find the last occurrence of 1
    last_one_index = len(input_grid) - 1 - input_grid[::-1].index(1)

    # Shift the sequence of 1s one position to the left
    output_grid = input_grid[:first_one_index] + input_grid[first_one_index+1:last_one_index+1] + [0] + input_grid[last_one_index+1:]
    return output_grid

# Test input
input_grid = [8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print("<<<" + " ".join(map(str, output_grid)) + ">>>")