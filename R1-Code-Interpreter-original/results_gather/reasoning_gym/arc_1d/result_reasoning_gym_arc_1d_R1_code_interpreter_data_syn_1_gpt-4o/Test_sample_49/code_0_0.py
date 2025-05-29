def transform_grid(input_grid):
    non_zero_elements = [x for x in input_grid if x != 0]
    zero_count = len(input_grid) - len(non_zero_elements)
    return non_zero_elements + [0] * zero_count

test_input = [0, 7, 7, 7, 7, 0, 0, 0, 0, 7, 7, 7, 7, 7, 0, 0, 0, 0, 7, 7, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)