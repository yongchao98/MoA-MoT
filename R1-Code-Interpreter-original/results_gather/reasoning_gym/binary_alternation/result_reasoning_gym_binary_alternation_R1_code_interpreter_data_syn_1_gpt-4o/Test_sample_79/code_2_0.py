def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Pattern 1: "010101..."
    swaps_0_start = 0
    swaps_1_start = 0
    
    for i in range(n):
        expected_char_0_start = '0' if i % 2 == 0 else '1'
        expected_char_1_start = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_0_start:
            swaps_0_start += 1
        if s[i] != expected_char_1_start:
            swaps_1_start += 1
    
    # If the string length is even, both patterns are possible
    if n % 2 == 0:
        return min(swaps_0_start // 2, swaps_1_start // 2)
    else:
        # If the string length is odd, only one pattern is possible
        if count_0 > count_1:
            return swaps_0_start // 2
        else:
            return swaps_1_start // 2

# Given binary string
binary_string = "100111101010101001"
result = min_swaps_to_alternate(binary_string)
print(result)