def find_max_iterative(nums):
    if len(nums) == 0:
        raise ValueError("find_max_iterative() arg is an empty sequence")
    max_num = nums[0]
    for x in nums:
        if x > max_num:
            max_num = x
    return {"max_value": max_num}

# Given input
input_data = {'nums': [9, 93, 85, -74, -84, -7, 78, 18, -5, -89, 55, 70, -57, 51]}
# Find the maximum value
result = find_max_iterative(input_data['nums'])
print(result)