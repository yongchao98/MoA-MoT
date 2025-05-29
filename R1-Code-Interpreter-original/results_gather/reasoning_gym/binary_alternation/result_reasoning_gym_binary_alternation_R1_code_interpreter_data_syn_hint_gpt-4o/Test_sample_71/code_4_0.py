def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed to match a given pattern
    def calculate_swaps(pattern):
        swaps = 0
        for i in range(n):
            if s[i] != pattern[i]:
                swaps += 1
        return swaps // 2
    
    # Generate the two possible patterns
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Calculate swaps for both patterns
    if n % 2 == 0:
        # If length is even, both patterns are possible
        swaps1 = calculate_swaps(pattern1)
        swaps2 = calculate_swaps(pattern2)
        return min(swaps1, swaps2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            return calculate_swaps(pattern1)
        else:
            return calculate_swaps(pattern2)

# Given binary string
binary_string = "10101110001100100100010101"
result = min_swaps_to_alternate(binary_string)
print(result)