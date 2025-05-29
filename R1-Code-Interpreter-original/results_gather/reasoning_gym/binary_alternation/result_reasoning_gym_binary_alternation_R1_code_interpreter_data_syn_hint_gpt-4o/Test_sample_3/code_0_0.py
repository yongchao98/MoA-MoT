def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to rearrange into an alternating pattern
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for both patterns
    mismatches_A = mismatches_B = 0
    for i in range(n):
        if i % 2 == 0:  # Even index
            if s[i] != '0':
                mismatches_A += 1
            if s[i] != '1':
                mismatches_B += 1
        else:  # Odd index
            if s[i] != '1':
                mismatches_A += 1
            if s[i] != '0':
                mismatches_B += 1
    
    # Calculate minimum swaps needed
    if n % 2 == 0:
        # If length is even, both patterns are possible
        min_swaps = min(mismatches_A, mismatches_B) // 2
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            min_swaps = mismatches_A // 2
        else:
            min_swaps = mismatches_B // 2
    
    return min_swaps

# Given binary string
binary_string = "101000001110001011001011111"
result = min_swaps_to_alternate(binary_string)
print(result)