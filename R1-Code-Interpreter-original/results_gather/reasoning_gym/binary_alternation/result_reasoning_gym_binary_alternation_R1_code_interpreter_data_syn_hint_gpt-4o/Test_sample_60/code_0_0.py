def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for pattern starting with '0'
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
    
    # Calculate minimum swaps
    if n % 2 == 0:
        # If length is even, both patterns are possible
        min_swaps = min(swaps_starting_with_0, swaps_starting_with_1) // 2
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            min_swaps = swaps_starting_with_0 // 2
        else:
            min_swaps = swaps_starting_with_1 // 2
    
    return min_swaps

# Given binary string
binary_string = "00101110001001111000110"
result = min_swaps_to_alternate(binary_string)
print(result)