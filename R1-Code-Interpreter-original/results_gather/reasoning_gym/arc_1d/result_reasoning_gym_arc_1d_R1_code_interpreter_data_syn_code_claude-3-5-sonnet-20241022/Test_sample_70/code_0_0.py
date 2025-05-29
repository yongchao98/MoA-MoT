def analyze_shift(input_str, output_str):
    # Convert strings to lists of integers
    input_list = [int(x) for x in input_str.split()]
    output_list = [int(x) for x in output_str.split()]
    
    # Find the different number and its positions
    unique_nums = set(input_list)
    shifts = {}
    
    for num in unique_nums:
        input_positions = [i for i, x in enumerate(input_list) if x == num]
        output_positions = [i for i, x in enumerate(output_list) if x == num]
        shift = output_positions[0] - input_positions[0]
        shifts[num] = shift
    
    return shifts

# Example inputs
example1_in = "0 0 0 0 0 9 9 9 9 9 9 9 9 9 9 9 9 0 0 0 0 0 0 0 0 0 0 0 0"
example1_out = "0 0 0 0 0 0 0 0 9 9 9 9 9 9 9 9 9 9 9 9 0 0 0 0 0 0 0 0 0"

example2_in = "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 0 0 0 0 0 0 3 3 3 3 3 3 3 3"
example2_out = "3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 0 0 0 0 0 0 3 3 3 3 3"

example3_in = "0 0 0 0 0 0 0 0 0 7 7 7 7 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 0"
example3_out = "0 0 0 0 0 0 0 0 0 0 0 0 7 7 7 7 7 7 7 7 7 7 7 7 7 0 0 0 0"

# Analyze shifts
print("Example 1 shifts:", analyze_shift(example1_in, example1_out))
print("Example 2 shifts:", analyze_shift(example2_in, example2_out))
print("Example 3 shifts:", analyze_shift(example3_in, example3_out))

# Test input
test_input = "1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1"

def apply_rule(input_str):
    nums = [int(x) for x in input_str.split()]
    n = len(nums)
    
    # Find the middle sequence (different number)
    main_num = nums[0]
    diff_num = None
    for num in nums:
        if num != main_num:
            diff_num = num
            break
    
    # Count sequences
    diff_count = nums.count(diff_num)
    
    # Create output
    result = [main_num] * n
    start_pos = 8  # Based on pattern analysis
    for i in range(diff_count):
        result[start_pos + i] = diff_num
    
    return ' '.join(map(str, result))

print("\nPredicted output:", apply_rule(test_input))