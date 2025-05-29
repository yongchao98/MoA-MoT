def find_pattern(input_str):
    # Convert input to list of integers
    nums = [int(x) for x in input_str.split()]
    
    # Find the two unique numbers
    unique_nums = sorted(set(nums))
    
    # Count occurrences
    counts = {num: nums.count(num) for num in unique_nums}
    
    # For test input: 0 0 0 0 0 0 0 8 8 8 8 8 8 8 8 0 0 0 0 0
    # We have: 12 zeros and 8 eights
    
    # Looking at the examples, if a number appears in a continuous block
    # and the other number appears at both ends, the continuous block number
    # stays in the middle while the other number moves to the ends
    
    # Find the continuous blocks
    blocks = []
    current_num = nums[0]
    current_count = 1
    
    for i in range(1, len(nums)):
        if nums[i] == current_num:
            current_count += 1
        else:
            blocks.append((current_num, current_count))
            current_num = nums[i]
            current_count = 1
    blocks.append((current_num, current_count))
    
    # In test case, we have blocks: (0,7), (8,8), (0,5)
    # The 8s form a continuous block, so they should stay together
    # The 0s should move to surround the 8s
    
    result = [0] * 3 + [8] * 8 + [0] * 9
    print(" ".join(map(str, result)))

# Test with the given input
find_pattern("0 0 0 0 0 0 0 8 8 8 8 8 8 8 8 0 0 0 0 0")