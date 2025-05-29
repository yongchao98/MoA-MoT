def analyze_pattern(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Count non-zeros
    input_nonzeros = sum(1 for x in input_list if x != 0)
    output_nonzeros = sum(1 for x in output_list if x != 0)
    
    # Find start and end positions of non-zeros
    input_start = next(i for i, x in enumerate(input_list) if x != 0)
    input_end = len(input_list) - next(i for i, x in enumerate(reversed(input_list)) if x != 0) - 1
    output_start = next(i for i, x in enumerate(output_list) if x != 0)
    output_end = len(output_list) - next(i for i, x in enumerate(reversed(output_list)) if x != 0) - 1
    
    print(f"Non-zeros count: Input={input_nonzeros}, Output={output_nonzeros}")
    print(f"Non-zero span: Input={input_start}-{input_end}, Output={output_start}-{output_end}")

# Analyze examples
print("Example 1:")
analyze_pattern("0 0 0 0 0 3 3 3 3 3 3 3 3 3 3 0", "0 3 3 3 3 3 3 3 3 3 3 0 0 0 0 0")
print("\nExample 2:")
analyze_pattern("0 0 0 0 3 3 3 3 3 3 3 3 3 3 0 0", "3 3 3 3 3 3 3 3 3 3 0 0 0 0 0 0")
print("\nExample 3:")
analyze_pattern("0 0 0 0 0 0 0 0 3 3 3 3 3 3 0 0", "0 0 0 0 3 3 3 3 3 3 0 0 0 0 0 0")