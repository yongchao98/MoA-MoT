def transform_grid(input_grid):
    # Find the start and end of the non-zero sequence
    start = next((i for i, x in enumerate(input_grid) if x != 0), None)
    end = next((i for i, x in enumerate(input_grid[::-1]) if x != 0), None)
    end = len(input_grid) - end - 1

    # Determine the shift based on the examples
    # From the examples, the shift is 2 positions to the left
    shift = 2

    # Create the output grid
    output_grid = [0] * len(input_grid)
    output_grid[start-shift:end-shift+1] = input_grid[start:end+1]

    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 9, 9, 9, 9]
test_output = transform_grid(test_input)

print("<<<" + " ".join(map(str, test_output)) + ">>>")