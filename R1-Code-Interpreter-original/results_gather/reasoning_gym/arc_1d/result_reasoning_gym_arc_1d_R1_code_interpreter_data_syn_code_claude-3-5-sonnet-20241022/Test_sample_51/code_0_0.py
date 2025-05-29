def analyze_pattern(input_str, output_str):
    input_nums = [int(x) for x in input_str.split()]
    output_nums = [int(x) for x in output_str.split()]
    
    # Function to find the rule
    def find_transformation_rule(inp, out):
        changes = []
        for i in range(len(inp)):
            if inp[i] != out[i]:
                changes.append((i, inp[i], out[i]))
        return changes

    # Analyze examples
    print("Transformation analysis:")
    changes = find_transformation_rule(input_nums, output_nums)
    for i, before, after in changes:
        print(f"Position {i}: {before} -> {after}")
        
    # Print consecutive sequences
    print("\nConsecutive 2s in input:")
    count = 0
    start = -1
    for i in range(len(input_nums)):
        if input_nums[i] == 2:
            if count == 0:
                start = i
            count += 1
        else:
            if count > 0:
                print(f"Sequence from {start} to {start+count-1}, length: {count}")
            count = 0
    if count > 0:
        print(f"Sequence from {start} to {start+count-1}, length: {count}")

# Analyze examples
print("Example 1:")
analyze_pattern("2 2 2 2 2 2 2 2 0 2 2 2 2 2 2 2 2 2 2 2 2 0 0 0",
                "8 8 8 8 2 2 2 2 0 8 8 8 8 8 8 2 2 2 2 2 2 0 0 0")

print("\nExample 2:")
analyze_pattern("2 2 0 0 2 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0",
                "8 2 0 0 8 8 8 8 8 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0")

print("\nExample 3:")
analyze_pattern("0 2 2 0 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0",
                "0 8 2 0 8 8 8 8 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0")