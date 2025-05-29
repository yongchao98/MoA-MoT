def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to form an alternating string
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Pattern 1: "010101..."
    mismatches_01 = sum(1 for i in range(n) if s[i] != str(i % 2))
    
    # Pattern 2: "101010..."
    mismatches_10 = sum(1 for i in range(n) if s[i] != str((i + 1) % 2))
    
    # Calculate swaps needed
    if n % 2 == 0:
        # If length is even, both patterns are possible
        swaps_01 = mismatches_01 // 2
        swaps_10 = mismatches_10 // 2
        return min(swaps_01, swaps_10)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            return mismatches_01 // 2
        else:
            return mismatches_10 // 2

# Given binary string
binary_string = "111101011010000001011"
result = min_swaps_to_alternate(binary_string)
print(result)