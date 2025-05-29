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
        if i % 2 == 0:  # Even index
            if s[i] != '0':
                swaps_starting_with_0 += 1
            if s[i] != '1':
                swaps_starting_with_1 += 1
        else:  # Odd index
            if s[i] != '1':
                swaps_starting_with_0 += 1
            if s[i] != '0':
                swaps_starting_with_1 += 1
    
    # If the string length is even, both patterns are possible
    # If the string length is odd, only one pattern is possible
    if n % 2 == 0:
        return min(swaps_starting_with_0, swaps_starting_with_1) // 2
    else:
        if count_0 > count_1:
            return swaps_starting_with_0 // 2
        else:
            return swaps_starting_with_1 // 2

# Given binary string
binary_string = "1100000110111"
result = min_swaps_to_alternate(binary_string)
print(result)