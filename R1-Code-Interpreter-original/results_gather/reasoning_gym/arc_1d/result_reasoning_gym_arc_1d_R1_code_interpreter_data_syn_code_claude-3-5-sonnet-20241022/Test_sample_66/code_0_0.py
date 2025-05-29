def analyze_pattern(input_grid, output_grid):
    # Convert string to list of integers
    input_nums = [int(x) for x in input_grid.split()]
    output_nums = [int(x) for x in output_grid.split()]
    
    # Find position of 2
    two_pos_in = input_nums.index(2)
    two_pos_out = output_nums.index(2)
    
    # Find sequence of non-zero numbers after 2
    sequence = []
    current_num = None
    count = 0
    
    for num in input_nums[two_pos_in+1:]:
        if num != 0 and num != 2:
            if current_num is None or num == current_num:
                current_num = num
                count += 1
            else:
                if count > 0:
                    sequence.append((current_num, count))
                current_num = num
                count = 1
    if count > 0:
        sequence.append((current_num, count))
    
    print(f"Position of 2: {two_pos_in}")
    print(f"Sequence after 2: {sequence}")
    print("Output arrangement:", output_nums)

# Test with Example 3
example3_in = "0 0 2 0 0 0 0 0 0 0 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8"
example3_out = "0 0 2 0 0 0 0 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 0 0 0"
print("Example 3 analysis:")
analyze_pattern(example3_in, example3_out)

# Test with our input
test_input = "0 2 0 0 0 0 0 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0"
print("\nTest input analysis:")
analyze_pattern(test_input, test_input)