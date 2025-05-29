def apply_final_rule(input_str):
    nums = input_str.split()
    i = 0
    while i < len(nums):
        if nums[i] == '0':
            # Found first zero, remove it
            nums.pop(i)
            nums.append('0')
            break
        i += 1
    return ' '.join(nums)

test_input = "3 8 4 0 0 0 7 7 1 0 5 0 0"
result = apply_final_rule(test_input)
print(result)