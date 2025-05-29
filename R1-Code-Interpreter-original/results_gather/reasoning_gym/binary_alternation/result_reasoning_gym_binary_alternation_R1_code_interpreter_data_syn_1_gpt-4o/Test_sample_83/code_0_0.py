def min_swaps_to_alternate(s):
    n = len(s)
    
    # Pattern 1: "010101..."
    mismatches_pattern1 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            mismatches_pattern1 += 1
    
    # Pattern 2: "101010..."
    mismatches_pattern2 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            mismatches_pattern2 += 1
    
    # Calculate swaps needed
    swaps_pattern1 = mismatches_pattern1 // 2 if mismatches_pattern1 % 2 == 0 else float('inf')
    swaps_pattern2 = mismatches_pattern2 // 2 if mismatches_pattern2 % 2 == 0 else float('inf')
    
    # Determine the minimum swaps
    min_swaps = min(swaps_pattern1, swaps_pattern2)
    
    # If min_swaps is inf, it means it's impossible to make the string alternating
    return min_swaps if min_swaps != float('inf') else -1

# Given binary string
binary_string = "1100100010111100101"
result = min_swaps_to_alternate(binary_string)
print(result)