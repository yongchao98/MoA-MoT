def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if first_non_zero is None:
        return input_grid  # If all are zeros, return as is

    # Find the last non-zero element
    last_non_zero = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1

    # Determine the block to shift
    block = input_grid[first_non_zero:last_non_zero + 1]

    # Shift the block 4 positions to the left
    new_position = max(0, first_non_zero - 4)
    output_grid = [0] * len(input_grid)
    output_grid[new_position:new_position + len(block)] = block

    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)