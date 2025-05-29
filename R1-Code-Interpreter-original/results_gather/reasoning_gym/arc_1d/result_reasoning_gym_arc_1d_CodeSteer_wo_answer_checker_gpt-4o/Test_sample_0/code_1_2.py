def find_shift(input_grid, output_grid):
    input_start = next((i for i, x in enumerate(input_grid) if x != 0), None)
    output_start = next((i for i, x in enumerate(output_grid) if x != 0), None)
    return input_start - output_start

def apply_rule(input_grid):
    # Based on the examples, the shift is consistently 2 positions to the left
    shift = 2
    non_zero_sequence = [x for x in input_grid if x != 0]
    output_grid = [0] * len(input_grid)
    output_grid[shift:shift+len(non_zero_sequence)] = non_zero_sequence
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 9, 9, 9, 9]
test_output = apply_rule(test_input)

print("<<<" + " ".join(map(str, test_output)) + ">>>")