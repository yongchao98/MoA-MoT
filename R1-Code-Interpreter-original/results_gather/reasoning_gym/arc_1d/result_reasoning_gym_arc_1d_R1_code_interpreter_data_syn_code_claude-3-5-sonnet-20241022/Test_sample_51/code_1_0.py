def transform_grid(input_str):
    nums = [int(x) for x in input_str.split()]
    result = nums.copy()
    
    # Process consecutive sequences
    count = 0
    start = -1
    for i in range(len(nums)):
        if nums[i] == 2:
            if count == 0:
                start = i
            count += 1
        else:
            if count > 0:
                # Apply the refined rule
                if count <= 2:
                    result[start] = 8
                elif count <= 8:
                    for j in range(start, start + 4):
                        result[j] = 8
                else:  # count >= 9
                    for j in range(start, start + 4):
                        result[j] = 8
            count = 0
            start = -1
    
    # Handle last sequence if it exists
    if count > 0:
        if count <= 2:
            result[start] = 8
        elif count <= 8:
            for j in range(start, start + 4):
                result[j] = 8
        else:  # count >= 9
            for j in range(start, start + 4):
                result[j] = 8
    
    return ' '.join(map(str, result))

test_input = "0 2 2 2 2 2 2 2 2 2 0 0 0 0 0 0 0 2 2 2 2 2 2 0"
output = transform_grid(test_input)
print(output)