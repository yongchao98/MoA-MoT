def transform_grid(input_grid):
    # Shift the grid left by one position
    shifted_grid = input_grid[1:]
    # Append a 0 to the end
    output_grid = shifted_grid + [0]
    return output_grid

test_input = [0, 3, 7, 9, 4, 0, 5, 0, 0, 0, 0, 2, 0, 3, 2, 3, 0]
test_output = transform_grid(test_input)
print("<<<" + " ".join(map(str, test_output)) + ">>>")