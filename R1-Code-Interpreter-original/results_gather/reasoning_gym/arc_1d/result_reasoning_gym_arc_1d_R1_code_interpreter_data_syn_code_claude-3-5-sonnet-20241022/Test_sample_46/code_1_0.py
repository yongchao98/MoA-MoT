def analyze_full_pattern(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Analyze non-zero numbers and their positions
    input_nonzeros = [(i, x) for i, x in enumerate(input_list) if x != 0]
    output_nonzeros = [(i, x) for i, x in enumerate(output_list) if x != 0]
    
    print("Non-zero numbers and their positions:")
    print(f"Input: {input_nonzeros}")
    print(f"Output: {output_nonzeros}")

print("Example 1:")
analyze_full_pattern("5 5 5 0 0 0 0 0 0 0 0 5 5", "5 5 0 0 0 0 0 0 0 0 5 5 5")
print("\nExample 2:")
analyze_full_pattern("1 1 1 1 1 1 1 0 1 1 1 1 1", "1 1 1 1 1 1 0 1 1 1 1 1 1")
print("\nExample 3:")
analyze_full_pattern("7 7 7 7 7 7 7 7 7 7 7 0 0", "7 7 7 7 7 7 7 7 7 7 0 0 7")