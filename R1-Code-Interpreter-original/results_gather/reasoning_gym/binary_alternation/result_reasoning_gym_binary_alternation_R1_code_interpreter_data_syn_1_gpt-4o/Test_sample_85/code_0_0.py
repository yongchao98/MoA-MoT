def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to rearrange into an alternating pattern
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
    
    # Each swap fixes two mismatches
    if count_0 == count_1:
        # If counts are equal, we can choose either pattern
        return min(swaps_starting_with_0, swaps_starting_with_1) // 2
    elif count_0 > count_1:
        # If there are more '0's, we must start with '0'
        return swaps_starting_with_0 // 2
    else:
        # If there are more '1's, we must start with '1'
        return swaps_starting_with_1 // 2

# Given binary string
binary_string = "11101011000101001"
result = min_swaps_to_alternate(binary_string)
print(result)