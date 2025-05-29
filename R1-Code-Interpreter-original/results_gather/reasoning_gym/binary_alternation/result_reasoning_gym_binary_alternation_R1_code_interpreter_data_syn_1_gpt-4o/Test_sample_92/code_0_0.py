def min_swaps_to_alternate(s):
    n = len(s)
    
    # Pattern 1: "010101..."
    mismatches1 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            mismatches1 += 1
    
    # Pattern 2: "101010..."
    mismatches2 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            mismatches2 += 1
    
    # Calculate swaps needed
    swaps1 = mismatches1 // 2 if mismatches1 % 2 == 0 else float('inf')
    swaps2 = mismatches2 // 2 if mismatches2 % 2 == 0 else float('inf')
    
    # Determine the minimum swaps
    min_swaps = min(swaps1, swaps2)
    
    return min_swaps if min_swaps != float('inf') else -1

# Given binary string
binary_string = "100001101101001011"
result = min_swaps_to_alternate(binary_string)
print(result)