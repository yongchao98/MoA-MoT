def analyze_transitions(input_grid, output_grid):
    # Analyze the pattern of transitions
    input_groups = []
    current_group = []
    current_num = None
    
    # Group consecutive non-zero numbers in input
    for num in input_grid:
        if num != 0 and (current_num is None or current_num == num):
            current_group.append(num)
            current_num = num
        elif num == 0 and current_group:
            input_groups.append(current_group)
            current_group = []
            current_num = None
    if current_group:
        input_groups.append(current_group)
    
    # Do the same for output
    output_groups = []
    current_group = []
    current_num = None
    for num in output_grid:
        if num != 0 and (current_num is None or current_num == num):
            current_group.append(num)
            current_num = num
        elif num == 0 and current_group:
            output_groups.append(current_group)
            current_group = []
            current_num = None
    if current_group:
        output_groups.append(current_group)
    
    print("Input groups:", input_groups)
    print("Output groups:", output_groups)

# Test with example 1
print("Example 1:")
example1_in = [0,6,0,0,6,6,6,6,6,6,0,0,6,6,6,0,0,6,6,0,6,6,6,0]
example1_out = [6,0,6,6,0,6,6,6,0,6,6,6,0,6,6,6,6,6,6,0,0,0,0,0]
analyze_transitions(example1_in, example1_out)