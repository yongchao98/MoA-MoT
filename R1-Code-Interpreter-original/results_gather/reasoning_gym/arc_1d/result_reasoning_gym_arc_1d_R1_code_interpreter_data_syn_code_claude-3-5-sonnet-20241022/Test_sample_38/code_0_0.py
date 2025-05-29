def find_pattern(input_str, examples):
    # Convert input strings to lists
    test_input = [int(x) for x in input_str.split()]
    
    # Find first non-zero number and count leading zeros
    non_zero_start = 0
    while non_zero_start < len(test_input) and test_input[non_zero_start] == 0:
        non_zero_start += 1
    
    # Find last non-zero number
    non_zero_end = len(test_input) - 1
    while non_zero_end >= 0 and test_input[non_zero_end] == 0:
        non_zero_end -= 1
    
    # Extract non-zero sequence
    non_zero_sequence = test_input[non_zero_start:non_zero_end + 1]
    
    # Create output with appropriate number of leading zeros
    output = [0] * non_zero_start + non_zero_sequence + [0] * (len(test_input) - len(non_zero_sequence) - non_zero_start)
    
    # Print result
    print(' '.join(map(str, output)))

# Test input
test_input = "0 0 0 0 6 3 5 1 7 9 1 2 4 9 6 3 2 0 0"
find_pattern(test_input, None)