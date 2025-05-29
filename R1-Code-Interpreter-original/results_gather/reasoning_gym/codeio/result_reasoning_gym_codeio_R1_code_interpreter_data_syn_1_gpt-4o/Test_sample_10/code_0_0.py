def find_max_iterative(nums):
    if len(nums) == 0:
        raise ValueError("find_max_iterative() arg is an empty sequence")
    max_num = nums[0]
    for x in nums:
        if x > max_num:
            max_num = x
    return {"max_value": max_num}

# Given input
input_data = {'nums': [-94, 89, -30, -38]}
# Find the maximum value
result = find_max_iterative(input_data['nums'])
print(result)