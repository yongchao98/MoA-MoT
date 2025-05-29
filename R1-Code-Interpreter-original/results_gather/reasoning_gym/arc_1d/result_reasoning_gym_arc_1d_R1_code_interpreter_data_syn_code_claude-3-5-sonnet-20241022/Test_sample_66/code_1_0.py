def analyze_redistribution(input_str, output_str):
    input_grid = [int(x) for x in input_str.split()]
    output_grid = [int(x) for x in output_str.split()]
    
    # Find the non-zero sequence
    non_zero = [(i, x) for i, x in enumerate(input_grid) if x != 0 and x != 2]
    sequence_num = non_zero[0][1] if non_zero else 0
    sequence_count = len(non_zero)
    
    # Find position of 2
    two_pos = input_grid.index(2)
    
    # Count leading and trailing zeros in output
    leading_zeros = 0
    trailing_zeros = 0
    for x in output_grid:
        if x == 0:
            leading_zeros += 1
        else:
            break
    for x in reversed(output_grid):
        if x == 0:
            trailing_zeros += 1
        else:
            break
            
    print(f"Number in sequence: {sequence_num}")
    print(f"Count in sequence: {sequence_count}")
    print(f"Position of 2: {two_pos}")
    print(f"Leading zeros in output: {leading_zeros}")
    print(f"Trailing zeros in output: {trailing_zeros}")
    print("Output pattern:", output_grid)

# Test with Example 3
example3_in = "0 0 2 0 0 0 0 0 0 0 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8"
example3_out = "0 0 2 0 0 0 0 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 0 0 0"
print("Example 3 analysis:")
analyze_redistribution(example3_in, example3_out)

# Test with our test input
test_input = "0 2 0 0 0 0 0 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0"
print("\nTest input analysis:")
analyze_redistribution(test_input, test_input)