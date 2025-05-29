def min_swaps_to_alternate(s):
    n = len(s)
    pattern1_mismatches = 0
    pattern2_mismatches = 0
    
    # Check mismatches for both patterns
    for i in range(n):
        expected_char_pattern1 = '0' if i % 2 == 0 else '1'
        expected_char_pattern2 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_pattern1:
            pattern1_mismatches += 1
        if s[i] != expected_char_pattern2:
            pattern2_mismatches += 1
    
    # Calculate swaps needed
    if pattern1_mismatches % 2 == 0:
        swaps_pattern1 = pattern1_mismatches // 2
    else:
        swaps_pattern1 = float('inf')  # Impossible to fix
    
    if pattern2_mismatches % 2 == 0:
        swaps_pattern2 = pattern2_mismatches // 2
    else:
        swaps_pattern2 = float('inf')  # Impossible to fix
    
    # Determine the minimum swaps needed
    min_swaps = min(swaps_pattern1, swaps_pattern2)
    
    return min_swaps if min_swaps != float('inf') else -1

# Given binary string
binary_string = "10101001100101101111100000"
result = min_swaps_to_alternate(binary_string)
print(result)