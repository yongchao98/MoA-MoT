def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to rearrange into an alternating pattern
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
    
    # Return the minimum swaps needed
    if count_0 == count_1:
        return min(swaps_start_with_0, swaps_start_with_1) // 2
    elif count_0 > count_1:
        return swaps_start_with_0 // 2
    else:
        return swaps_start_with_1 // 2

# Given binary string
binary_string = "001101001110100110010101011100"
result = min_swaps_to_alternate(binary_string)
print(result)