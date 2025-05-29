def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed for a given pattern
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
    swaps1 = calculate_swaps(pattern1)
    swaps2 = calculate_swaps(pattern2)
    
    # If counts of 0s and 1s are equal, return the minimum of both swaps
    if count_0 == count_1:
        return min(swaps1, swaps2)
    # If there are more 0s, pattern1 is the only valid pattern
    elif count_0 > count_1:
        return swaps1
    # If there are more 1s, pattern2 is the only valid pattern
    else:
        return swaps2

# Given binary string
binary_string = "111110000100"
# Calculate the minimum swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)