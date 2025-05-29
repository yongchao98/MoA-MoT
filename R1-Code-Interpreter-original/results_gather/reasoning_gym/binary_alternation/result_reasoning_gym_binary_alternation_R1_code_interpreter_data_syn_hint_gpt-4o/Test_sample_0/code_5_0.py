def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for pattern starting with '0' (e.g., "010101...")
    swaps_starting_with_0 = 0
    swaps_starting_with_1 = 0
    
    for i in range(n):
        expected_char_0 = '0' if i % 2 == 0 else '1'
        expected_char_1 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_0:
            swaps_starting_with_0 += 1
        if s[i] != expected_char_1:
            swaps_starting_with_1 += 1
    
    # We need to divide by 2 because each swap fixes two positions
    return min(swaps_starting_with_0 // 2, swaps_starting_with_1 // 2)

# Given binary string
binary_string = "01000111101110"
result = min_swaps_to_alternate(binary_string)
print(result)