def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference in count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Pattern 1: "010101..."
    mismatches_01 = sum(1 for i in range(n) if s[i] != str(i % 2))
    
    # Pattern 2: "101010..."
    mismatches_10 = sum(1 for i in range(n) if s[i] != str((i + 1) % 2))
    
    # Calculate swaps needed for each pattern
    swaps_01 = mismatches_01 // 2
    swaps_10 = mismatches_10 // 2
    
    # If the string length is even, both patterns are possible
    if n % 2 == 0:
        return min(swaps_01, swaps_10)
    else:
        # If the string length is odd, only one pattern is possible
        if count_0 > count_1:
            return swaps_01
        else:
            return swaps_10

# Given binary string
binary_string = "110000101101101110001001"
result = min_swaps_to_alternate(binary_string)
print(result)