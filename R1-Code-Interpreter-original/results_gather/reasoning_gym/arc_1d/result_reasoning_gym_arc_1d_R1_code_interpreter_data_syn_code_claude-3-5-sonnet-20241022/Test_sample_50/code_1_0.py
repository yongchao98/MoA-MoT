def analyze_groups(input_str, output_str):
    input_arr = [int(x) for x in input_str.split()]
    output_arr = [int(x) for x in output_str.split()]
    
    # Find groups in input
    input_groups = []
    current_group = []
    for i, x in enumerate(input_arr):
        if x != 0:
            current_group.append(i)
        elif current_group:
            input_groups.append(current_group)
            current_group = []
    if current_group:
        input_groups.append(current_group)
    
    # Find groups in output
    output_groups = []
    current_group = []
    for i, x in enumerate(output_arr):
        if x != 0:
            current_group.append(i)
        elif current_group:
            output_groups.append(current_group)
            current_group = []
    if current_group:
        output_groups.append(current_group)
    
    print("Input groups:", [len(g) for g in input_groups])
    print("Input group positions:", input_groups)
    print("Output groups:", [len(g) for g in output_groups])
    print("Output group positions:", output_groups)
    print("---")

# Test input
test_input = "2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2 2"
print("Test case:")
analyze_groups(test_input, "")