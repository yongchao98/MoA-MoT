def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for pattern starting with '0' and '1'
    swaps_start_with_0 = swaps_start_with_1 = 0
    
    # Pattern starting with '0': "010101..."
    for i in range(n):
        if i % 2 == 0 and s[i] != '0':
            swaps_start_with_0 += 1
        elif i % 2 == 1 and s[i] != '1':
            swaps_start_with_0 += 1
    
    # Pattern starting with '1': "101010..."
    for i in range(n):
        if i % 2 == 0 and s[i] != '1':
            swaps_start_with_1 += 1
        elif i % 2 == 1 and s[i] != '0':
            swaps_start_with_1 += 1
    
    # If the number of 0s and 1s are equal, we can choose either pattern
    if count_0 == count_1:
        return min(swaps_start_with_0, swaps_start_with_1) // 2
    # If there are more 0s, we must start with '0'
    elif count_0 > count_1:
        return swaps_start_with_0 // 2
    # If there are more 1s, we must start with '1'
    else:
        return swaps_start_with_1 // 2

# Given binary string
binary_string = "011001010101010011101001010110"
result = min_swaps_to_alternate(binary_string)
print(result)