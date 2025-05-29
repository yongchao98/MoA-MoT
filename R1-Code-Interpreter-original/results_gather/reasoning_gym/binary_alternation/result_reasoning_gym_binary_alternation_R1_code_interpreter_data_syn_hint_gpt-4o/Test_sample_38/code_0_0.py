def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check feasibility
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to count swaps needed for a given pattern
    def count_swaps(pattern):
        swaps = 0
        for i in range(n):
            if s[i] != pattern[i]:
                swaps += 1
        return swaps // 2
    
    # Generate possible patterns
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Calculate swaps for each pattern
    if n % 2 == 0:
        # If length is even, both patterns are possible
        swaps1 = count_swaps(pattern1)
        swaps2 = count_swaps(pattern2)
        return min(swaps1, swaps2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            return count_swaps(pattern1)
        else:
            return count_swaps(pattern2)

# Given binary string
binary_string = "1100000110111"
result = min_swaps_to_alternate(binary_string)
print(result)