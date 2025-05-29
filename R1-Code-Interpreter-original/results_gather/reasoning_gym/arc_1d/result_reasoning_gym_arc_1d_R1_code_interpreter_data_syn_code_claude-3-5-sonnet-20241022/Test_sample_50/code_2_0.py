def analyze_transformation(input_str, output_str):
    input_arr = [int(x) for x in input_str.split()]
    output_arr = [int(x) for x in output_str.split()]
    
    # Find groups and their sizes
    def find_groups(arr):
        groups = []
        current_size = 0
        for x in arr:
            if x != 0:
                current_size += 1
            elif current_size > 0:
                groups.append(current_size)
                current_size = 0
        if current_size > 0:
            groups.append(current_size)
        return groups
    
    input_groups = find_groups(input_arr)
    output_groups = find_groups(output_arr)
    
    print(f"Input groups sizes: {input_groups}")
    print(f"Output groups sizes: {output_groups}")
    print(f"Total numbers: {sum(input_groups)} -> {sum(output_groups)}")
    print("---")

# Analyze all examples
print("Example 1:")
analyze_transformation("3 0 0 0 0 0 0 0 0 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3",
                      "3 3 3 3 0 0 0 0 0 0 0 0 3 3 3 3 3 3 3 3 3 3 3 3 3")

print("Example 2:")
analyze_transformation("5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 5 5",
                      "5 5 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5")

print("Example 3:")
analyze_transformation("1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1",
                      "0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")

print("Test case:")
analyze_transformation("2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2", "")