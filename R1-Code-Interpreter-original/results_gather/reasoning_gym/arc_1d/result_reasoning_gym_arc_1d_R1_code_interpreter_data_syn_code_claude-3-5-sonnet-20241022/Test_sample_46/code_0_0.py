def analyze_pattern(input_str, output_str):
    # Convert strings to lists
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find positions of zeros
    input_zeros = [i for i, x in enumerate(input_list) if x == 0]
    output_zeros = [i for i, x in enumerate(output_list) if x == 0]
    
    print(f"Input zeros at: {input_zeros}")
    print(f"Output zeros at: {output_zeros}")
    
    # Print the shifts
    if input_zeros and output_zeros:
        print(f"Shift: {[out - inp for inp, out in zip(input_zeros, output_zeros)]}")

# Test examples
print("Example 1:")
analyze_pattern("5 5 5 0 0 0 0 0 0 0 0 5 5", "5 5 0 0 0 0 0 0 0 0 5 5 5")
print("\nExample 2:")
analyze_pattern("1 1 1 1 1 1 1 0 1 1 1 1 1", "1 1 1 1 1 1 0 1 1 1 1 1 1")
print("\nExample 3:")
analyze_pattern("7 7 7 7 7 7 7 7 7 7 7 0 0", "7 7 7 7 7 7 7 7 7 7 0 0 7")