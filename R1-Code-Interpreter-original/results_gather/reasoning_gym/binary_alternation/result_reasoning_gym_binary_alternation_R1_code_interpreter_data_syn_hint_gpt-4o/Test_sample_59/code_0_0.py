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
    
    # Check feasibility
    if pattern1_mismatches % 2 != 0 and pattern2_mismatches % 2 != 0:
        return -1
    
    # Calculate minimum swaps
    swaps_pattern1 = pattern1_mismatches // 2
    swaps_pattern2 = pattern2_mismatches // 2
    
    # Return the minimum swaps
    return min(swaps_pattern1, swaps_pattern2)

# Given binary string
binary_string = "101111000100"
# Calculate the minimum swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)