def min_swaps_to_alternate(s):
    n = len(s)
    
    # Pattern 1: Starts with '0'
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    # Pattern 2: Starts with '1'
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Count mismatches for both patterns
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    # Calculate swaps needed
    swaps1 = mismatches1 // 2 if mismatches1 % 2 == 0 else float('inf')
    swaps2 = mismatches2 // 2 if mismatches2 % 2 == 0 else float('inf')
    
    # Determine the minimum swaps
    min_swaps = min(swaps1, swaps2)
    
    return min_swaps if min_swaps != float('inf') else -1

# Given binary string
binary_string = "101001001010000111011"
result = min_swaps_to_alternate(binary_string)
print(result)