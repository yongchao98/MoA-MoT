def min_swaps_to_alternate(s):
    n = len(s)
    pattern1 = "01" * (n // 2) + "0" * (n % 2)
    pattern2 = "10" * (n // 2) + "1" * (n % 2)
    
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    # Calculate swaps needed
    swaps1 = mismatches1 // 2 if mismatches1 % 2 == 0 else float('inf')
    swaps2 = mismatches2 // 2 if mismatches2 % 2 == 0 else float('inf')
    
    # Determine the minimum swaps needed
    min_swaps = min(swaps1, swaps2)
    
    return min_swaps if min_swaps != float('inf') else -1

# Given binary string
binary_string = "010001001011111110100101"
print(min_swaps_to_alternate(binary_string))