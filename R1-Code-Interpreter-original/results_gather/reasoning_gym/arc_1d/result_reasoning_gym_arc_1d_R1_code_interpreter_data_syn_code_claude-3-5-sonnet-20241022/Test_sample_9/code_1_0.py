def analyze_sequences(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find first and last non-zero numbers
    first_nonzero_in = next((i for i, x in enumerate(input_list) if x != 0), -1)
    last_nonzero_in = len(input_list) - 1 - next((i for i, x in enumerate(reversed(input_list)) if x != 0), -1)
    
    # Find first and last non-zero numbers in output
    first_nonzero_out = next((i for i, x in enumerate(output_list) if x != 0), -1)
    last_nonzero_out = len(output_list) - 1 - next((i for i, x in enumerate(reversed(output_list)) if x != 0), -1)
    
    print(f"First non-zero: {first_nonzero_in} -> {first_nonzero_out}")
    print(f"Last non-zero: {last_nonzero_in} -> {last_nonzero_out}")
    
    # Count changes
    changes = sum(1 for i in range(len(input_list)) if input_list[i] != output_list[i])
    print(f"Total changes: {changes}")

print("Example 1:")
analyze_sequences("0 0 7 7 7 7 7 7 7 7 7 0 7 0 0 7 0 0 7 0", 
                 "0 0 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 0 0 0")

print("\nExample 2:")
analyze_sequences("1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1", 
                 "0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1")

print("\nExample 3:")
analyze_sequences("0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0", 
                 "0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0")