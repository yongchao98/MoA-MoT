def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check feasibility
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps for both patterns
    swaps_start_with_0 = 0
    swaps_start_with_1 = 0
    
    for i in range(n):
        expected_char_0 = '0' if i % 2 == 0 else '1'
        expected_char_1 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_0:
            swaps_start_with_0 += 1
        if s[i] != expected_char_1:
            swaps_start_with_1 += 1
    
    # If the string length is even, both patterns are possible
    if n % 2 == 0:
        return min(swaps_start_with_0, swaps_start_with_1) // 2
    else:
        # If the string length is odd, only one pattern is possible
        if count_0 > count_1:
            return swaps_start_with_0 // 2
        else:
            return swaps_start_with_1 // 2

# Given binary string
binary_string = "11100100001101110"
result = min_swaps_to_alternate(binary_string)
print(result)