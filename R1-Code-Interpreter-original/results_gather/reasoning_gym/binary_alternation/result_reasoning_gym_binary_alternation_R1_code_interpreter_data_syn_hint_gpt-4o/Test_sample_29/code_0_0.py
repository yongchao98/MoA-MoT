def min_swaps_to_alternate(s):
    n = len(s)
    pattern_a = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern_b = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    mismatches_a = sum(1 for i in range(n) if s[i] != pattern_a[i])
    mismatches_b = sum(1 for i in range(n) if s[i] != pattern_b[i])
    
    # Calculate the minimum swaps needed
    swaps_a = mismatches_a // 2
    swaps_b = mismatches_b // 2
    
    # Check if it's possible to make the string alternating
    if mismatches_a % 2 != 0 and mismatches_b % 2 != 0:
        return -1
    
    return min(swaps_a, swaps_b)

# Given binary string
binary_string = "1100101011000010110101011"
result = min_swaps_to_alternate(binary_string)
print(result)