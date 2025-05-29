def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for pattern "010101..."
    swaps_0_start = 0
    swaps_1_start = 0
    
    for i in range(n):
        if i % 2 == 0:  # Even index
            if s[i] != '0':
                swaps_0_start += 1
            if s[i] != '1':
                swaps_1_start += 1
        else:  # Odd index
            if s[i] != '1':
                swaps_0_start += 1
            if s[i] != '0':
                swaps_1_start += 1
    
    # If the string length is even, both patterns are possible
    # If the string length is odd, only one pattern is possible
    if n % 2 == 0:
        return min(swaps_0_start, swaps_1_start) // 2
    else:
        if count_0 > count_1:
            return swaps_0_start // 2
        else:
            return swaps_1_start // 2

# Given binary string
binary_string = "0100110100110100101101001110"
result = min_swaps_to_alternate(binary_string)
print(result)