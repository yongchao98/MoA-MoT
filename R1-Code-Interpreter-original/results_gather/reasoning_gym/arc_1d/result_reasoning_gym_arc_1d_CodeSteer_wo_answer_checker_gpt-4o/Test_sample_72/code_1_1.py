def shift_sequence(input_grid, shift_by=3):
    # Find the start and end of the non-zero sequence
    start, end = None, None
    for i, num in enumerate(input_grid):
        if num != 0:
            if start is None:
                start = i
            end = i

    # If no non-zero sequence is found, return the original grid
    if start is None or end is None:
        return input_grid

    # Calculate the new start position after shifting
    new_start = start + shift_by

    # Create a new grid filled with zeros
    output_grid = [0] * len(input_grid)

    # Place the non-zero sequence in the new position
    for i in range(start, end + 1):
        if new_start < len(output_grid):
            output_grid[new_start] = input_grid[i]
            new_start += 1

    return output_grid

# Test the function with the given input
test_input = [0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = shift_sequence(test_input)
print("<<<" + " ".join(map(str, output)) + ">>>")