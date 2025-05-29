def analyze_pattern(input_str, output_str):
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Print differences
    print("Differences at positions:")
    for i in range(len(input_list)):
        if input_list[i] != output_list[i]:
            print(f"Position {i}: {input_list[i]} -> {output_list[i]}")
    
    # Analyze pattern
    changes = []
    for i in range(len(input_list)):
        if input_list[i] != output_list[i]:
            changes.append((i, input_list[i], output_list[i]))
    
    return changes

# Example 1
ex1_in = "0 0 7 7 7 7 7 7 7 7 7 0 7 0 0 7 0 0 7 0"
ex1_out = "0 0 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 0 0 0"
print("Example 1 changes:")
analyze_pattern(ex1_in, ex1_out)

# Example 2
ex2_in = "1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1"
ex2_out = "0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1"
print("\nExample 2 changes:")
analyze_pattern(ex2_in, ex2_out)

# Example 3
ex3_in = "0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 0"
ex3_out = "0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0"
print("\nExample 3 changes:")
analyze_pattern(ex3_in, ex3_out)