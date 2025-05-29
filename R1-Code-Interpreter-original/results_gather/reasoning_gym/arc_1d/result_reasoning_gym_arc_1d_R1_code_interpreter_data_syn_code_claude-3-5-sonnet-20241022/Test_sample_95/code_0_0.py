def shift_pattern(input_grid):
    # Convert input string to list of integers
    nums = [int(x) for x in input_grid.split()]
    
    # Find first non-zero number
    start = 0
    while start < len(nums) and nums[start] == 0:
        start += 1
    
    # Find last non-zero number
    end = len(nums) - 1
    while end >= 0 and nums[end] == 0:
        end -= 1
        
    # If all zeros or no zeros, return original
    if start > end:
        return nums
    
    # Create result array with zeros
    result = [0] * len(nums)
    
    # Shift non-zero sequence right by 3
    for i in range(start, end + 1):
        result[i + 3] = nums[i]
    
    return result

# Test input
test = "0 0 0 0 0 0 0 0 0 0 0 0 0 3 2 8 9 4 2 5 4 7 7 2 9 7 0 0 0"

# Get result and format output
result = shift_pattern(test)
print(" ".join(map(str, result)))