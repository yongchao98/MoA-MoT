def apply_correct_rule(input_str):
    nums = input_str.split()
    i = 0
    while i < len(nums) - 1:
        if nums[i] == '0':
            # If next number is also 0, remove one 0
            if nums[i + 1] == '0':
                nums.pop(i)
                nums.append('0')
                i += 1  # Skip the next 0
            # If it's a single 0, remove it
            else:
                nums.pop(i)
                nums.append('0')
        else:
            i += 1
    return ' '.join(nums)

test_input = "3 8 4 0 0 0 7 7 1 0 5 0 0"
result = apply_correct_rule(test_input)
print(result)