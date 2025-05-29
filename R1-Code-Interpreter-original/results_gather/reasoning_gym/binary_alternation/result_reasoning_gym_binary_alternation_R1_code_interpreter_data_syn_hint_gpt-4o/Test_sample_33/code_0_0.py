def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for pattern "010101..."
    mismatches_pattern_a = sum(1 for i in range(n) if s[i] != str(i % 2))
    
    # Calculate mismatches for pattern "101010..."
    mismatches_pattern_b = sum(1 for i in range(n) if s[i] != str((i + 1) % 2))
    
    # If the string length is even, both patterns are possible
    if n % 2 == 0:
        return min(mismatches_pattern_a, mismatches_pattern_b) // 2
    else:
        # If the string length is odd, only one pattern is possible
        if count_0 > count_1:
            return mismatches_pattern_a // 2
        else:
            return mismatches_pattern_b // 2

# Given binary string
binary_string = "00111110010111001000001"
result = min_swaps_to_alternate(binary_string)
print(result)