def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps for pattern starting with '0'
    swaps_0 = 0
    swaps_1 = 0
    for i in range(n):
        if i % 2 == 0:  # Even index
            if s[i] != '0':
                swaps_0 += 1
            if s[i] != '1':
                swaps_1 += 1
        else:  # Odd index
            if s[i] != '1':
                swaps_0 += 1
            if s[i] != '0':
                swaps_1 += 1
    
    # Choose the minimum swaps needed
    if n % 2 == 0:
        # If length is even, both patterns are possible
        return min(swaps_0, swaps_1) // 2
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            return swaps_0 // 2
        else:
            return swaps_1 // 2

# Given binary string
binary_string = "00111100111010101111000000"
result = min_swaps_to_alternate(binary_string)
print(result)