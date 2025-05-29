def analyze_changes(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find first non-zero position in input and output
    first_nonzero_in = next((i for i, x in enumerate(input_list) if x != 0), len(input_list))
    first_nonzero_out = next((i for i, x in enumerate(output_list) if x != 0), len(output_list))
    
    # Find last non-zero position in input and output
    last_nonzero_in = len(input_list) - 1 - next((i for i, x in enumerate(reversed(input_list)) if x != 0), len(input_list))
    last_nonzero_out = len(output_list) - 1 - next((i for i, x in enumerate(reversed(output_list)) if x != 0), len(output_list))
    
    print(f"First non-zero shift: {first_nonzero_out - first_nonzero_in}")
    print(f"Last non-zero shift: {last_nonzero_out - last_nonzero_in}")
    print(f"Non-zero length change: {(last_nonzero_out - first_nonzero_out) - (last_nonzero_in - first_nonzero_in)}")

print("Example 1:")
analyze_changes("0 0 0 0 0 0 0 6 6 6 6 6 6 6 6 6 6 6 6 6 6", 
               "0 0 0 0 0 0 6 6 6 6 6 6 6 6 6 6 6 6 6 6 0")
print("\nExample 2:")
analyze_changes("0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0",
               "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0")
print("\nExample 3:")
analyze_changes("0 0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0",
               "0 0 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0 0")