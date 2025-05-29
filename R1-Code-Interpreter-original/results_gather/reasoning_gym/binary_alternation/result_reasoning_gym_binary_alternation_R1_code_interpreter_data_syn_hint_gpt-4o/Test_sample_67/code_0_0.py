def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to rearrange into an alternating pattern
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps for pattern starting with '0' (i.e., "010101...")
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
    
    # Calculate minimum swaps needed
    if count_0 == count_1:
        # If counts are equal, we can choose either pattern
        min_swaps = min(swaps_starting_with_0, swaps_starting_with_1) // 2
    elif count_0 > count_1:
        # More '0's, must start with '0'
        min_swaps = swaps_starting_with_0 // 2
    else:
        # More '1's, must start with '1'
        min_swaps = swaps_starting_with_1 // 2
    
    return min_swaps

# Given binary string
binary_string = "1001000101010101110110"
result = min_swaps_to_alternate(binary_string)
print(result)