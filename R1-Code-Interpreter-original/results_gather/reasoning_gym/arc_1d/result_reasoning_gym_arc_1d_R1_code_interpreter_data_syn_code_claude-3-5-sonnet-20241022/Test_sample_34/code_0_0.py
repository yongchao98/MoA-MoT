def analyze_pattern(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find non-zero sequences
    def get_sequence(lst):
        start = next((i for i, x in enumerate(lst) if x != 0), None)
        if start is None:
            return None, None, None, 0
        end = next((i for i, x in enumerate(lst[start:], start) if x == 0), len(lst))
        return start, end, lst[start], end - start

    in_start, in_end, value, length = get_sequence(input_list)
    out_start, out_end, _, _ = get_sequence(output_list)
    
    shift = in_start - out_start
    
    # For test input
    test_input = "0 0 0 0 0 0 9 9 9 9 9 9 0 0 0 0 0 0 0"
    test_list = [int(x) for x in test_input.split()]
    test_start, test_end, test_value, test_length = get_sequence(test_list)
    
    # Apply the derived rule
    result = [0] * 19
    new_start = test_start - 2  # Shift left by 2 positions
    for i in range(test_length):
        result[new_start + i] = test_value
    
    print(" ".join(map(str, result)))

# Test with Example 1
input1 = "0 0 0 0 8 8 8 8 8 8 8 8 8 8 8 8 8 8 0"
output1 = "0 0 8 8 8 8 8 8 8 8 8 8 8 8 8 8 0 0 0"
analyze_pattern(input1, output1)