def min_swaps_to_alternate(s):
    n = len(s)
    pattern_a = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern_b = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    mismatches_a = sum(1 for i in range(n) if s[i] != pattern_a[i])
    mismatches_b = sum(1 for i in range(n) if s[i] != pattern_b[i])
    
    return min(mismatches_a, mismatches_b)

binary_string = "01001010110101010"
result = min_swaps_to_alternate(binary_string)
print(result)