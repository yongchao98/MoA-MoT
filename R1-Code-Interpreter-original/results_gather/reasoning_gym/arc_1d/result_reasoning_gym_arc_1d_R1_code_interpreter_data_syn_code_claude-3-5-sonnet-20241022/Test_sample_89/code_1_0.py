def find_rule(input_str):
    input_list = [int(x) for x in input_str.split()]
    
    # Find non-zero sequence properties
    nonzero_count = sum(1 for x in input_list if x != 0)
    first_nonzero = next(i for i, x in enumerate(input_list) if x != 0)
    last_nonzero = len(input_list) - next(i for i, x in enumerate(reversed(input_list)) if x != 0) - 1
    
    # Calculate output positions
    output_start = max(0, 16 - (last_nonzero - first_nonzero + 1) - 5)
    output_end = output_start + nonzero_count - 1
    
    # Generate output
    output = ['0'] * 16
    non_zero_value = input_list[first_nonzero]
    for i in range(output_start, output_start + nonzero_count):
        output[i] = str(non_zero_value)
    
    print("Input analysis:")
    print(f"Non-zero count: {nonzero_count}")
    print(f"Input span: {first_nonzero}-{last_nonzero}")
    print(f"Output span: {output_start}-{output_end}")
    print("Generated output:", ' '.join(output))

# Test with the given input
test_input = "0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 0"
find_rule(test_input)