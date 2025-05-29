def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to rearrange into an alternating pattern
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for both patterns
    swaps_start_with_0 = swaps_start_with_1 = 0
    
    # Pattern starting with '0'
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            swaps_start_with_0 += 1
    
    # Pattern starting with '1'
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            swaps_start_with_1 += 1
    
    # Since each swap fixes two positions, divide by 2
    if count_0 == count_1:
        return min(swaps_start_with_0, swaps_start_with_1) // 2
    elif count_0 > count_1:
        return swaps_start_with_0 // 2
    else:
        return swaps_start_with_1 // 2

# Given binary string
binary_string = "1111111000010000"
result = min_swaps_to_alternate(binary_string)
print(result)