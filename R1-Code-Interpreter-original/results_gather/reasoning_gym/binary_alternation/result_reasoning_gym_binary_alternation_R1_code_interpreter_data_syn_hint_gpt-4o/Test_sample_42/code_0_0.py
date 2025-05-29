def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps for pattern starting with '0'
    swaps_starting_with_0 = 0
    swaps_starting_with_1 = 0
    
    for i in range(n):
        expected_char_0 = '0' if i % 2 == 0 else '1'
        expected_char_1 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_0:
            swaps_starting_with_0 += 1
        if s[i] != expected_char_1:
            swaps_starting_with_1 += 1
    
    # Since each swap fixes two positions, we divide by 2
    if count_0 == count_1:
        return min(swaps_starting_with_0, swaps_starting_with_1) // 2
    elif count_0 > count_1:
        return swaps_starting_with_0 // 2
    else:
        return swaps_starting_with_1 // 2

# Given binary string
binary_string = "00100110101100110010010111101"
result = min_swaps_to_alternate(binary_string)
print(result)