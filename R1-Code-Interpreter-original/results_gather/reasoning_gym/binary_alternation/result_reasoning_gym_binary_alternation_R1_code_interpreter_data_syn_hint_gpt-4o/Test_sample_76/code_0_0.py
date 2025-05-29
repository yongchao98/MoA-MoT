def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed for a given pattern
    def calculate_swaps(s, pattern):
        swaps = 0
        for i in range(len(s)):
            if s[i] != pattern[i]:
                swaps += 1
        return swaps // 2
    
    # Generate the two possible patterns
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Calculate swaps for both patterns
    swaps1 = calculate_swaps(s, pattern1)
    swaps2 = calculate_swaps(s, pattern2)
    
    # Determine the minimum swaps
    if count_0 == count_1:
        return min(swaps1, swaps2)
    elif count_0 > count_1:
        return swaps1
    else:
        return swaps2

# Given binary string
binary_string = "10001010101101"
result = min_swaps_to_alternate(binary_string)
print(result)